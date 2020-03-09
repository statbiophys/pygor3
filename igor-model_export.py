#!/home/alfaceor/anaconda3/envs/pygor3/bin/python

import pygor3 as p3
import pandas as pd

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

# task.run_evaluate()
task.update_batch_filenames()
task.load_IgorModel()

for strEvent in task.mdl.xdata.keys():
    da = task.mdl.xdata[strEvent]
    dependencias = list(task.mdl.xdata[strEvent].dims)
    dependencias.remove(strEvent)
    dependencias_dim = [task.mdl.xdata[strEvent][dep].shape[0] for dep in dependencias]
    print(strEvent, da.dims, dependencias_dim, dependencias)
    if len(dependencias) == 0 :
        lbl_file = "P__"+strEvent
        print(lbl_file)
    elif len(dependencias) == 1:
        lbl_file = "P__" + strEvent +"__g__"+dependencias[0]
        print(lbl_file)
        df = pd.DataFrame(data=da.values, index=da['lbl__' + dependencias[0]].values,
                          columns=da['lbl__' + strEvent].values)
        #print(df)
    elif len(dependencias) == 2 :
        lbl_file = "P__" + strEvent + "__g__" + dependencias[0]+"__"+ dependencias[1]
        print(lbl_file)
    elif len(dependencias) == 3:
        lbl_file = "P__" + strEvent + "__g__" + dependencias[0] + "__" + dependencias[1]+ "__" + dependencias[2]
        print(lbl_file)
    else:
        print("OHHHHHHHH")
    #df = pd.DataFrame(data=da.values, index=da['lbl__' + 'v_choice'].values, columns=da['lbl__' + 'j_choice'].values)


task.run_clean_batch()
