#!/usr/bin/env python3
#
# Kim Brugger (17.08.2022) kbr(at)brugger.dk

import sys
import os
import pprint as pp
import fcsparser
import flowkit as fk
from bokeh.plotting import show
import matplotlib.pyplot as plt
_ = plt.ioff()

def rev_dict(input:dict) -> dict:
    res = {}
    for k in input.keys():
        res[ input[k]] = k

    return res


def pns_to_pnn(pns:list, pnn:list) -> dict:
    res = {}
    for i in range(0, len(pns)):
        res[pns[i]] = pnn[i]

    return res

def main():

    filename = 'processed/barcode_1/fcs/A2_specific.fcs'
#    filename = 'raw/barcode_12/Barcode_12_12_Processed.fcs'
 #   meta = fcsparser.parse(filename, meta_data_only=True)
 #   print( type(meta))
 #   print( meta.keys())
 
 #   sys.exit()

#    meta, data = fcsparser.parse(filename, meta_data_only=False, reformat_meta=True)
 #   print( type(meta))
#    print( meta.keys())
#    print(meta['__header__'])
 #   print(meta['_channels_'])

#    print( data )

    sample = fk.Sample(filename, subsample=10000)
    print( sample )
    pp.pprint(sample.channels )
    meta = sample.get_metadata()
    rev_meta = pns_to_pnn(sample.pns_labels, sample.pnn_labels)
    pp.pprint(rev_meta)
#    sys.exit()
#    pp.pprint( meta)
#    print( meta['Pt198Di'])
#    print( "PNN:", sample.pnn_labels )
#    print( "PNS:", sample.pns_labels )
#    print( sample.get_channel_index('CD4'))
#    print( "Rev cd16: ", rev_meta['CD16'])
#    print( sample.get_channel_number_by_label('CD16'))
#    print( "Rev 195Pt: ", rev_meta['195Pt']) 
#    print(sample.get_channel_number_by_label('Pt195Di'))
#    print(sample.get_channel_index('195Pt'))
#    print( sample.get_channel_number_by_label(rev_meta['195Pt']))
#    print( "Subsample size: ", len(sample.subsample_indices) )

#    pp.pprint( rev_meta )

#    p = sample.plot_histogram('Bi209Di', source='raw', bins=256)
#    show(p)

#    f = sample.plot_channel('Bi209Di', source='raw')
#    plt.show()
#    sys.exit()
#    f = sample.plot_contour('Ir191Di', 'Ir193Di', source='raw')
#    plt.show()

#  yttrium (89Y), indium (115In),
# cerium (140Ce), terbium (159Tb), lutetium (175Lu), and bismuth (209Bi).

# 'p20s': '115In',
# 'p27s': '140Ce',
# 'p46s': '159Tb_CD45RO',
# 'p62s': '175Lu_CD10',
# 'p73s': '209Bi_CD16',
# 'p8s': '89Y_CD45',



    xform = fk.transforms.LogicleTransform('asinh', param_t=1024, param_w=0.5, param_m=4.5, param_a=0)
    sample.apply_transform(xform)

    f = sample.plot_contour(rev_meta['CD20'], rev_meta['CD19'], source='xform', plot_events=True)
    plt.show()





    for bead_channel in ['115In','140Ce','159Tb_CD45RO','175Lu_CD10','209Bi_CD16','89Y_CD45',]:
        try:
            print('plotting ', bead_channel)
            f = sample.plot_contour('Time', rev_meta[bead_channel], source='raw', plot_events=True)
            plt.show()
        except Exception as e:
            print( 'failed: ', e)


#    f = sample.plot_contour(rev_meta['CD16'], rev_meta['CD56'], source='raw', plot_events=True)
#    plt.show()

#    sys.exit()

#    CD16, CD56, CD3, CD38

#    cd4 = sample.get_channel_events(sample.get_channel_index(rev_meta['CD3']), source='raw')
#    cd3 = sample.get_channel_events(sample.get_channel_index(rev_meta['CD14']), source='raw')
#    print( cd4 )
#    for i, _ in enumerate(cd4):
#        print(cd3[i], cd4[i])

#    sample.as_dataframe(source='xform')

    channels = ['Time', 'Event_length', 'Center', 'Width', 'Residual', 'Offset', 'Amplitude', ]
    channels_index = []
    for channel in channels:
        channels_index.append( sample.get_channel_index(channel) )

    print("\t".join(map(str, channels)))

    for event in sample.get_events(source='raw', subsample=False):
        obs = []
        for channel_index in channels_index:
            obs.append( event[ channel_index])
#            print( obs )
        print("\t".join(map(str, obs)))



 #   cd3_index = sample.get_channel_index(rev_meta['CD3'])
 #   cd16_index = sample.get_channel_index(rev_meta['CD16'])
 #   cd38_index = sample.get_channel_index(rev_meta['CD38'])
 #   cd56_index = sample.get_channel_index(rev_meta['CD56'])

 #   killers = 0
 #   total   = 0
 #   for event in sample.get_events(source='raw', subsample=False):
 #       total += 1
 #       if min(event[cd3_index], event[cd16_index],event[cd38_index],event[cd56_index]) > 0.01:
 #           killers += 1

#    print( f"{killers} killer-cells in {total} observations, {killers/total}")



#get_channel_number_by_label
if __name__ == "__main__":
    main()
