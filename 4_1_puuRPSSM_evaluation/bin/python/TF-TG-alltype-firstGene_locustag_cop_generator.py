from collections import defaultdict
import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog= 'TF-TG-alltype-firstGene_locustag_cop_generator', 
                                     description="This script finds the head of the operon of a set of TGs", 
                                     epilog='For suggestions/questions please refer to ericka.hernandez@ibt.unam.mx')
    parser.add_argument('-p','--operon_path', type=str, help='File path where the operon prediction is')
    parser.add_argument('-i','--interactions_path', type=str, help='File with the TF-TG interactions is')
    parser.add_argument('-o','--output_path', type=str, help='Outputh path where the final files should be saved')
    args = parser.parse_args()

# operon_path = "/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/docs/eco.ope"
# interactions_path = "/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/docs/TF-TG-alltype-firstGene_locustag.tsv"
# output_path = "/space24/PGC/emhernan/4_1_puuRPSSM_evaluation/docs/TF-TG-alltype-firstGene_locustag_cop.tsv"

    operons = {}
    table = defaultdict(list)
    i = 0
    
    with open(args.operon_path, 'r') as f:
        for line in f:
            genes = line.strip().split()
            if genes:
                cop = genes[0]
                for gen in genes:
                    operons[gen] = cop
    
    with open(args.interactions_path, 'r')  as f:
        for line in f:
            if i != 0:
                fields = line.strip().split('\t')
                table['riType'].append(fields[0])
                table['tfName'].append(fields[1])
                table['tgName'].append(fields[2])
                table['riEvidence'].append(fields[3])
                table['tgLocusTag'].append(fields[4])
                key = 'eco-{}'.format(fields[4])
                if key in operons:
                    try:
                        table['tgCopLocusTag'].append(operons[key].split('-')[1])
                    except IndexError:
                        table['tgCopLocusTag'].append("NaN")
                else:
                    table['tgCopLocusTag'].append("NaN")
            i+=1
    
    df = pd.DataFrame(table)
    df_reordered = df[['riType', 'tfName', 'tgName', 'riEvidence', 'tgLocusTag','tgCopLocusTag']]
    df_reordered.to_csv(args.output_path, sep="\t", index = False)