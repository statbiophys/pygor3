#!/home/alfaceor/anaconda3/envs/pygor3/bin/python

import pygor3 as p3

task = p3.IgorTask()

from  optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--species", dest="species", help='Igor species')
parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

(options, args) = parser.parse_args()

if options.species == None:
    print("species is mandatory")
    exit()
else:
    task.igor_specie = str(options.species)

if options.chain == None:
    print("species is mandatory")
    exit()
else:
    task.igor_chain = str(options.chain)

if options.output == None:
    print("output is mandatory")
    exit()
else:
    output = options.output

if options.batch == None:
    print("Batchname not set. Random batchname will be assing.")

if len(args) > 0 and len(args) <2:
  filename = str(args[0])
  task.igor_read_seqs = filename
else:
  print("Supply a correct filename")
  exit()

task.run_evaluate()

# task = p3.IgorTask.load_from_batchname("testing")
# task.igor_specie = "human"
# task.igor_chain = "beta"
# task.load_IgorModel()
task.update_batch_filenames()
task.load_IgorModel()
df = task.get_pgen_pd()
task.run_clean_batch()
df.to_csv(output )