import pandas as pd

# load the GFF3 gene annotation and the predicted
annotations = pd.read_csv("C:\\Users\\tinas\\Downloads\\Vibrio_vulnificus.ASM74310v1.37.gff3", 
                          skiprows = 394, header = None, sep = '\t')
gff_file = pd.read_csv('results.gff3', header = None, sep = '\t')

# filter the GFF3 annotations to get CDS features on the forward strand
annotations = annotations[annotations.iloc[:, 2] == "CDS"]
annotations = annotations[annotations.iloc[:, 6] == "+"]

# store the starting and stopping positions for both the annotated and predicted 
annotated_starts = list(annotations.iloc[:, 3].values)
annotated_stops = list(annotations.iloc[:, 4].values)
predicted_starts = list(gff_file.iloc[:, 3].values)
predicted_stops = list(gff_file.iloc[:, 4].values)

annotated_match_both, annotated_match_start, annotated_match_end, annotated_no_match = 0, 0, 0, 0

# keep track of counts of the annotated matches
for i in range (0, len(annotated_starts), 1):
    if (int(annotated_starts[i]) in predicted_starts) and (int(annotated_stops[i]) in predicted_stops):
        annotated_match_both += 1
    elif int(annotated_starts[i]) in predicted_starts:
        annotated_match_start += 1
    elif int(annotated_stops[i]) in predicted_stops:
        annotated_match_end += 1
    else:
        annotated_no_match += 1

# get the fraction of the matches for annotated
annotated_match_both = float(annotated_match_both)/float((len(annotated_starts)))
annotated_match_start = float(annotated_match_start)/float((len(annotated_starts)))
annotated_match_end = float(annotated_match_end)/float((len(annotated_starts)))
annotated_no_match = float(annotated_no_match)/float((len(annotated_starts)))

print("Fractions based on annotated genes:")
print("Perfectly match both ends: " + str(annotated_match_both * 100))
print("Match the start only: " + str(annotated_match_start * 100))
print("Match the end only: " + str(annotated_match_end * 100))
print("Match neither: " + str(annotated_no_match * 100))

predicted_match_both, predicted_match_start, predicted_match_end, predicted_no_match = 0, 0, 0, 0

# keep track of counts of the predicted matches
for i in range (0, len(predicted_starts), 1):
    if (float(predicted_starts[i]) in annotated_starts) and (float(predicted_stops[i]) in annotated_stops) :
        predicted_match_both += 1
    elif float(predicted_starts[i]) in annotated_starts:
        predicted_match_start += 1
    elif float(predicted_stops[i]) in annotated_stops:
        predicted_match_end += 1
    else:
        predicted_no_match += 1

# get the fraction of the matches for predicted
predicted_match_both = float(predicted_match_both)/float((len(predicted_starts)))
predicted_match_start = float(predicted_match_start)/float((len(predicted_starts)))
predicted_match_end = float(predicted_match_end)/float((len(predicted_starts)))
predicted_no_match = float(predicted_no_match)/float((len(predicted_starts)))

print("Fractions based on predicted genes:")
print("Perfectly match both ends: " + str(predicted_match_both * 100))
print("Match the start only: " + str(predicted_match_start * 100))
print("Match the end only: " + str(predicted_match_end * 100))
print("Match neither: " + str(predicted_no_match * 100))
