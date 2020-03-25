#!/usr/bin/env python3
import pygor3 as p3
import pandas as pd

from optparse import OptionParser

def main():
    parser = OptionParser()

    (options, args) = parser.parse_args()
    # fln_model_parms = args[0]
    # fln_model_marginals = args[1]
    # mdl = p3.IgorModel(model_parms_file=fln_model_parms, model_marginals_file=fln_model_marginals)
    mdl = p3.IgorModel.load_default('human', 'tcr_beta')
    # for strEvent in mdl.xdata.keys():
    #     mdl.xdata[strEvent]
    strEvent = 'v_choice'
    print(mdl.xdata[strEvent].dims)
    print('*'*50)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    strEvent = 'v_choice'
    df = mdl.xdata[strEvent].to_dataframe(name='P')
    df.plot.bar(x='lbl__' + strEvent, y='P', rot=90, ax=ax)
    fig.tight_layout()
    fig.savefig(strEvent + ".pdf")

    fig, ax = plt.subplots()
    strEvent = 'j_choice'
    da = mdl.xdata[strEvent]
    da = da.sum(dim='v_choice')
    P_J = da
    df = da.to_dataframe(name='P')

    df.plot.bar(x='lbl__' + strEvent, y='P', rot=90, ax=ax)
    fig.tight_layout()
    fig.savefig(strEvent + ".pdf")

    fig, ax = plt.subplots()
    strEvent='d_gene'
    da = mdl.xdata[strEvent] * mdl.xdata['j_choice']
    da = da.sum(dim='j_choice')
    da.plot(ax=ax, marker='o', label=str(i))

#mdl_anchs = p3.IgorAnchors(flnVanchors, flnJanchors)


if __name__ == "__main__":
    main()
