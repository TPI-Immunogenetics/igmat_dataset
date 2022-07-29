import os
import re
import sys
import yaml
import time
import urllib.request
import urllib.parse
import urllib.error
import argparse
import traceback
from subprocess import Popen, PIPE

import src.imgt as imgt
import src.fasta as fasta

base_path  = os.path.split(__file__)[0]

def generateAlignment(data, species, path):

  # Delete file if exists
  if os.path.exists(path):
    os.remove(path)

  # # Initialize file
  # with open(path, 'w') as handle:
  #   handle.write("# STOCKHOLM 1.0\n")

  for chain in data:
    if 'V' not in data[chain] or 'J' not in data[chain]:
      raise Exception('Missing V or J region')

    headerSize = 0
    sequenceSize = 0
    result = []
    for v in data[chain]['V']:
      for j in data[chain]['J']:

        # Combine the sequence
        combined = v['sequence'] + j['sequence']

        # Store it
        name = '{species}|{v_allele}|{j_allele}'.format(
          species=species,
          v_allele=v['name'],
          j_allele=j['name']
        ).replace(" ", "_")

        headerSize = max(headerSize, len(name))
        sequenceSize = max(sequenceSize, len(combined))
        result.append({
          'name': name,
          'sequence': combined
        })

    with open(path, 'a') as handle:
      handle.write("# STOCKHOLM 1.0\n")
      handle.write("#=GF ID {species}_{chain}\n".format(
        species=species.replace(' ', '+'),
        chain=chain
      ))
      for i in range(len(result)):
        handle.write('{name} {sequence}\n'.format(
          name=result[i]['name'].replace(' ', '_').ljust(headerSize+1),
          sequence=result[i]['sequence'].replace('.', '-')
        ))

      handle.write("{end} {match}\n".format(end="#=GC RF".ljust(headerSize+1), match="x"*sequenceSize))
      handle.write("//\n")

def validateAlignment(path, verbose=False):

  ruleList = {
    23: {
      'residues': ['C'],
      'critic': True
    },
    41: {
      'residues': ['W'],
      'critic': False
    },
    89: {
      'residues': ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V', 'P'],
      'critic': False
    },
    104: {
      'residues': ['C'],
      'critic': True
    },
    118: {
      'residues': ['F', 'W'],
      'critic': False
    }
  }

  try:

    warningCount = 0
    if not os.path.exists(path):
      raise Exception('Unable to find alignment file')

    chainList = []
    filename = os.path.basename(path)
    with open(path, 'r') as handle:

      count = 0
      chain_type = None
      for line in handle:

        count += 1
        line = line.strip()

        # Parse ID line
        if line.startswith('#=GF ID'):
          chain_type = line.replace('#=GF ID', '').split('_')[1]
          chainList.append(chain_type)
          continue

        # Skip comment lines
        if not line or line[0] == '#':
          continue

        # End of alignment
        if line.startswith('//'):
          type = None
          continue

        # Extract sequence data
        line = line.split()
        if len(line) != 2:
          raise Exception('Invalid format at line %d' % count)

        head = line[0]
        sequence = line[1]
        for index in ruleList:

          try:
                    
            if sequence[index-1] not in ruleList[index]['residues']:

              message = 'Chain {0},line {1}: position {2} must be one of \'{3}\''.format(
                chain_type,
                count,
                index,
                ','.join(ruleList[index]['residues'])
                )
              raise (Exception(message) if ruleList[index]['critic'] else Warning(message))

          except Warning as w:

            warningCount += 1
            if verbose:
              print('Warning {0}: {1}'.format(filename, w))
              index_min = max(0, index-1-10)
              index_max = min(index-1+10, len(sequence))
              print('  ' + '.. ' + sequence[index_min:index_max] + ' ..')
              print('  ' + '   ' + (' ' * (index-1-index_min) + '*'))
              # print('  ' + sequence)
              # print('  ' + (' ' * (index-1) + '*'))
            

  except Exception as e:
    print('Error validating alignment {0}: {1}'.format(filename, str(e)))
    return False

  if warningCount:
    print('Found {0} warnings while validating data. Please check the alignment'.format(warningCount))

  return chainList
  # return True

def loadConfig(path):

  ''' Load config file and check for validity '''
  if not os.path.exists(path):
    raise Exception('Unable to find config file')

  with open(config_path, 'r') as handle:
    config = yaml.safe_load(handle)

  # Check config file
  if 'species' not in config:
    raise Exception('No species section')

  if 'chains' not in config:
    raise Exception('No chain url list')

  return config

def fetchAlignment(url, output_path, force=False, retry=3):
  '''
  Fetch data from IMGT and store it in cache
  '''
  
  # Check if the file already exists
  if force or not os.path.exists(output_path):

    # Download data from IMGT
    count = retry
    data = None
    while count > 0:
      try:
        handle = urllib.request.urlopen(url)
        data = handle.read()

      except Exception as e:
        # print('Unable to retrieve url %s [retry %d/%d]' % (url, count, retry))
        traceback.print_exc()
        count -= 1
        time.sleep(0.1)
        data = None
        continue

      # Data downloaded successfully
      break;

    # Something happens
    if data is None:
      return False

    # Parse extracted html data
    parser = imgt.IMGTDBParser()
    sequences = parser.rip_sequences(data.decode('utf-8'))
    if not sequences:
      return False
      # raise Exception('No sequences available for URL \'{url}\''.format(url=url))
    
    # All good, store fasta file
    with open(output_path, "w") as handle:
      for name, sequence in sequences:
        # Skip invalid sequences
        if not name:
          continue
          
        # Cleanup name
        name = name[1:] if name[0] == '>' else name

        # Store fasta
        handle.write('>{name}\n{sequence}\n'.format(name=name, sequence=sequence))

  imgt_fields = [
    "accession_number",
    "allele",  
    "species",  
    "functionality",  
    "region",  
    "start_and_end_positions_IMGT/LIGM-DB_accession_number", 
    "number_of_nucleotides", 
    "codon_start", 
    "number_nucleotides_added_in_5'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_added_in_3'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_to_correct_sequencing_errors", 
    "number_of_amino_acids", 
    "number_of_characters", 
    "partial",  
    "reverse"
  ]

  # Load sequence data
  result = []
  for sequence in fasta.parse(output_path):

    # Parse sequence name
    fields = [x.strip() for x in sequence.getName().split("|")]
    fields = dict(list(zip(imgt_fields, fields)))

    # Filter out sequences
    isPartial = True if len(fields['partial']) > 0 else False
    isReverse = True if len(fields['reverse']) > 0 else False
    isGermline = True if fields["allele"].split("*")[-1].strip() == "01" else False

    try:

      # Check sequence validity
      if not sequence.isValid(gap=True):
        raise Warning('Invalid sequence')

      # Skip sequence with no accession number
      if fields['accession_number'] == 'None':
        raise Warning('Unknown accession number')
      
      # Skip non functional, incomplete or reverse sequences
      if fields["functionality"] != "F" or isPartial or isReverse:
        raise Warning('Non-functional [functionality=\'{0}\'], partial [partial=\'{1}\'] or reverse [reverse=\'{2}\'] sequence'.format(
          fields['functionality'], 
          fields['partial'],
          fields['reverse']
        ))

      if not isGermline:
        raise Warning('Not a germline sequence: {0}'.format(fields['allele']))

      # All good, store sequence
      result.append({
        'name': fields['allele'],
        'sequence': sequence.getSequence()
      })
    except Warning as e:
      print('Skipping allele {0}: {1}'.format(fields['allele'], str(e)))

  return result

if __name__ == "__main__":

  # Initialize the argument parser
  argparser = argparse.ArgumentParser()
  argparser.add_argument('--force', '-f', action="store_true", dest="force", help="force rebuild")
  argparser.add_argument('--verbose', '-v', action="store_true", dest="verbose", help="Increase verbosity")
  
  # Parse arguments
  args = argparser.parse_args();

  # Check dist folder
  dist_path = os.path.join(base_path, 'dist')
  if not os.path.exists(dist_path):
    os.mkdir(dist_path)
    
  # Check cache folder
  cache_path = os.path.join(base_path, 'cache')
  if not os.path.exists(cache_path):
    os.mkdir(cache_path)

  # Load config file
  config_path = os.path.join(base_path, 'config.yml')
  config = loadConfig(config_path)

  # Load reference from IMGT
  resultList = []
  for species in config['species']:

    print('Processing species {0}'.format(species))
    output_filename = '{0}.sto'.format(species).replace(' ', '_').lower()
    output_stockholm = os.path.join(dist_path, output_filename)
    if not os.path.exists(output_stockholm) or args.force:
      
      data = {}
      for chain in config['chains']:

        # Prepare the url
        chain_type = chain[0]
        region_type = chain[1]
        name = species.replace(' ', '+')
        url = config['chains'][chain].format(species=name)
        try:
        
          # Extract data
          fasta_path = os.path.join(cache_path, '{0}_{1}.fasta'.format(species, chain).replace(' ', '_').lower())
          alignment = fetchAlignment(url, fasta_path)
          if not alignment:
            raise Exception('Unable to fetch {0} for organism {1}'.format(chain, species))

          # Store alignment data
          data.setdefault(chain_type, {})[region_type] = alignment

        except Exception as e:
          print('\tUnable to get {0} {1}: {2}'.format(species, chain, str(e)))

      # Format IMGT genes
      print('Formatting IMGT data')
      data = imgt.format(data, species)

      # Generate alignments
      generateAlignment(data, species, output_stockholm)
  
    # Check the alignment
    chainList = validateAlignment(output_stockholm, verbose=args.verbose)
    
    # Update result list
    resultList.append({
      'species': species,
      'path': output_filename,
      'valid': True if chainList else False,
      'chain': chainList
    })

  # Store results 
  list_path = os.path.join(dist_path, 'list.txt')
  with open(list_path, 'w') as handle:

    handle.write('species\tpath\tchain\n')
    for i in range(len(resultList)):
      if not resultList[i]['valid']:
        print('Alignment {0} is not valid. Please check'.format(resultList[i]['species']))
        continue

      handle.write('{0}\t{1}\t{2}\n'.format(
        resultList[i]['species'],
        resultList[i]['path'],
        ','.join(resultList[i]['chain'])
      ))



    
  
