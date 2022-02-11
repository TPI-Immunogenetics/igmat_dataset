
from html.parser import HTMLParser
from html.entities import name2codepoint

import os
import sys
import shutil
import traceback
import tempfile
from subprocess import Popen, PIPE

from . import fasta
# from . import helpers


file_path  = os.path.split(__file__)[0]
bin_path = os.path.abspath(os.path.join(file_path, "../bin"))

class IMGTDBParser(HTMLParser):
  currenttag = None
  currentnamedent = None
  _data = []

  def handle_starttag(self, tag, attrs):
    self.currenttag=tag

  def handle_endtag(self, tag):
    self.currenttag=None

  def handle_data(self, data):
      split = data.split("\n")
      start = sum([ 1 if l[0]==">" else 0 for l in split if len(l)])
      if self.currenttag=="pre" and (self.currentnamedent ==">" or start):
          # Two different ways of parsing the html based on how IMGT have formatted the pages.
          # For some reason they format gene db differently sometimes (legacy?) 
          if start > 1: # If you encounter more than one line in the data with a fasta ">" symbol, all sequences will be in the same packet
              name, sequence = None, ""
              for l in split:
                  if not l: continue
                  if l[0]==">":
                      if sequence:
                          self._data.append( (name, sequence) )
                          name, sequence = None, ""
                      name = l
                  else:
                      sequence += l.replace(" ", "")
              if name and sequence:
                  self._data.append( (name, sequence) )
          else: # Otherwise it will be done entry by entry
              print("1")
              try:
                  name = split[0]
              except IndexError:
                  return
              sequence = ("".join( split[1:])).replace(" ", "")
              self._data.append( (name, sequence) )

  def handle_entityref(self, name):
      self.currentnamedent = chr(name2codepoint[name])

  def handle_charref(self, name):
      if name.startswith('x'):
          self.currentnamedent = chr(int(name[1:], 16))
      else:
          self.currentnamedent = chr(int(name))

  def rip_sequences(self,htmlstring):
      """
      Method for this subclass that automates the return of data
      """
      self.reset()
      self._data = []
      self.currenttag = None
      self.currentnamedent = None
      self.feed(htmlstring)
      self._data
      return self._data

def _mouse_delta(sequence):
  """
  Mouse delta chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.

  This is particularly bad because alignment is not even consistent within the chain and species!!!

  Remove and return
  """
  # Check in here because not all are bad...recheck again in the format v genes just to make sure.
  if sequence[103] != "C" or sequence[22] != "C":
      return sequence[ : 8 ] + sequence[ 9:85 ] + sequence[86:]
  return sequence

def _rhesus_lambda(sequence):
  """
  Rhesus lambda chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
  Remove and return
  """
  return sequence[:20]+sequence[21:51]+ sequence[53:] 

def _mouse_alpha(sequence):
  """
  Mouse alpha chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
  Remove and return
  """
  return sequence[:8]+sequence[9:85]+sequence[86:]

def _format_j_genes(jalignments, species, chain):

  # The reference sequence 
  reference = {
    'name': 'Mus|H|Mus musculus|IGHJ3*01',
    'sequence': 'WFAYWGQGTLVTVSA',
    'start': 4,
    'stop': 19,
  }

  # Get alignment executable
  # musclePath = helpers.get_dir_binary('muscle')
  musclePath = 'muscle'

  # Process sequences
  # results = {}
  results = []
  tempFolder = tempfile.mkdtemp()
  try:

    # Write a fasta file containing all sequences
    hasReference = False
    input_file = os.path.join(tempFolder, 'input.fa')
    output_file = os.path.join(tempFolder, 'output.fa')
    with open(input_file, "w" ) as handle:
      # for al in jalignments:
      for i in range(len(jalignments)):

        # Trim down sequence to 15 characters
        trimmed = jalignments[i]['sequence'][-15:]
        handle.write('>{species}|{chain}|{type}|{allele}\n{sequence}\n'.format(
          species=species,
          chain=chain,
          type='V',
          allele=jalignments[i]['name'],
          sequence=trimmed
        ))

      # If the reference sequence is not in the list, add it
      if not hasReference:
        handle.write('>{name}\n{sequence}\n'.format(
          name=reference['name'],
          sequence=reference['sequence']
        ))

    # Run muscle
    process = Popen( [ musclePath, "-in", input_file, "-gapopen", "-10", "-out", output_file, ], stdout=PIPE, stderr=PIPE )
    _, pr_stderr = process.communicate()

    if not os.path.exists(output_file):
      raise Exception('Error while executing Muscle.')

    # Extract reference sequence
    for sequence in fasta.parse(output_file):
      if sequence.getName() == reference['name']:
        ref_aligned = sequence.getSequence()
        break

    # start = ref_aligned.index(reference['sequence'])
    # start = 0
    start = max(0, ref_aligned.find(reference['sequence']))
    START = (start+1-reference['start']) if start > reference['start'] else 0
    END = start + len(reference['sequence'])
    for sequence in fasta.parse(output_file):

      if not hasReference and sequence.getName() == reference['name']:
        continue

      species, chain, chain_type, allele = sequence.getName().strip(">").split("|")

      # We take the last 13 of the new alignment and pad into 20 long string 
      padded = sequence.getSequence()[START: END][-14:].rjust(20).replace(" ", ".")
      results.append({
        'name': allele,
        'sequence': padded
      })

  finally:

    # Remove temp folder
    shutil.rmtree(tempFolder)

  return results

def _format_v_genes(valignments, species, chain):

  # These are special functions for specific chains
  speciesVFormat = {
    'Macaca+mulatta_LV': _rhesus_lambda,
    'Mus_AV': _mouse_alpha,
    'Mus_DV': _mouse_delta
  }
  
  results = {} 
  for entry in valignments:
    results = []
    for data in valignments:

      try:

        # sequence = valignments[entry][seq]
        sequence = data['sequence']
        name = '{species}_{type}{chain}'.format(species=species, type=chain, chain='V')
        if name in speciesVFormat:
          sequence = speciesVFormat[name](sequence)

        # Trim sequence to the right size (108) and pad with gaps on the right side
        results.append({
          'name': data['name'],
          'sequence': sequence[:108].ljust(108).replace(" ",".")
        })

      except Exception as e:

        print('Unable to format V genes: %s' % e)

  return results

def format(alignments, species):

  # for chain in alignments:
  for chain in list(alignments.keys()):

    try:

      # Format V region
      if 'V' not in alignments[chain]:
        raise Exception('Unable to find V-regions for chain {0} in {1}'.format(chain, species))

      alignments[chain]['V'] = _format_v_genes(alignments[chain]['V'], species, chain)

      # Format J region
      if 'J' not in alignments[chain]:
        raise Exception('Unable to find J-regions for chain {0} in {1}'.format(chain, species))

      alignments[chain]['J'] = _format_j_genes(alignments[chain]['J'], species, chain)
    except Exception as e:
      print('Discarding chain {0} in {1}: {2}'.format(chain, species, str(e)))
      alignments.pop(chain)

  return alignments