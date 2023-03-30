
from Bio import SeqIO
from Bio import Seq

########################################################################
import sys
import argparse
# geneious IO shit
parser=argparse.ArgumentParser()

parser.add_argument("-input", help="input option")
parser.add_argument("-output", help="output option")

args=parser.parse_args()

print(f"Args: {args}\nCommand Line: {sys.argv}\ninput: {args.input}")
print(f"Dict format: {vars(args)}")


gen_obj = SeqIO.read(str(args.input), "genbank")

##########################################################################

class Guide:
    def __init__(self, label, start, end, strand, sequence):
        self.label = label
        self.start = start
        self.end = end
        self.strand = strand
        self.sequence = sequence
        self.len = self.end - self.start



# parse out info
gene = gen_obj.seq
translation = gene.translate(to_stop=True)
feature = gen_obj.features

# candidate sites for editing
pos_sites = ['CAA', 'CAG', 'CGA']
pos_coords = []
neg_sites = ['TGG']
neg_coords = []

window_dict = {18: (16,18),
               20: (17,19),
               22: (18,20)}


# iterate through the gene identifying coordinates that can be mutated to stop codons
for i in range(0, len(gene), 3):
    if gene[i:i+3] in pos_sites:
        pos_coords.append(i)
    elif gene[i:i+3] in neg_sites:
        neg_coords.append(i)
    else:
        pass

# scrape the list of guide sequences
pos_scraped = [Guide(str(feature.qualifiers['label'][0]), int(feature.location.start), int(feature.location.end), int(feature.strand), str(feature.qualifiers['Target_Sequence'][0])) for feature in gen_obj.features if feature.type == 'CRISPR' and feature.strand == 1]
neg_scraped = [Guide(str(feature.qualifiers['label'][0]), int(feature.location.start), int(feature.location.end), int(feature.strand), str(feature.qualifiers['Target_Sequence'][0])) for feature in gen_obj.features if feature.type == 'CRISPR' and feature.strand == -1]

# create secondary list for guides which overlap with the coordinates of interest
candidates = []
candidate_seqs = []


for guide in pos_scraped:
    window = range(((guide.end - 3) - window_dict[guide.len-3][1]), ((guide.end-3) - window_dict[guide.len-3][0]) + 1)

    for coord in pos_coords:
        if coord in window:
            candidates.append(guide)
            candidate_seqs.append(guide.sequence)
            break

for guide in neg_scraped:
    window = range(((guide.start + 3) + window_dict[guide.len - 3][0]) - 1, ((guide.start + 3) + window_dict[guide.len - 3][1]))
    for coord in neg_coords:
        if (coord + 1) in window and (coord + 2) in window:
            candidates.append(guide)
            candidate_seqs.append(guide.sequence)
            break

# save only the desired features, rewrite to the object that will be output
#candidate_features = [feat for feat in gen_obj.features if feat.type == 'source' or feat.type == 'CDS' or (feat.type == 'CRISPR' and feat.qualifiers['label'][0] in candidate_labels)]
candidate_features = [feat for feat in gen_obj.features if feat.type == 'source' or feat.type == 'CDS' or (feat.type == str('CRISPR') and str(feat.qualifiers['Target_Sequence'][0]) in candidate_seqs)]

#candidate_features = [feat for feat in gen_obj.features if (not (feat.type == 'CRISPR' and feat.qualifiers['Target_Sequence'][0] in candidate_seqs)) or (feat.type == 'CRISPR' and feat.qualifiers['Target_Sequence'][0] in candidate_seqs)]

gen_obj.features = candidate_features
############################################

# write to file
with open(str(args.output), "w") as output_handle:
    SeqIO.write(gen_obj, output_handle, "genbank")




