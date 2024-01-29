import sys
import csv

csv.field_size_limit(sys.maxsize)

class Item:
  def __init__(self, chr=0, src='', feature='', start=0, end=0, score=0, strand='', frame='', attribute=''):
    self.chr = chr
    self.src = src
    self.feature = feature
    self.start = start
    self.end = end 
    self.score = score
    self.strand = strand
    self.frame = frame
    self.attribute = attribute

class Group:
    def __init__(self, circRNA=0, exon_list=[]):
        self.circRNA = circRNA
        self.exon_list = exon_list

def input_to_groups(input):

    groups = []

    #read each line, form a group of circRNA with exons and add to groups
    with open(input, "r", encoding="utf8") as circ_star_file:
        tsv_reader = csv.reader(circ_star_file, delimiter="\t")

        #Skip the first row, which is the header
        #next(tsv_reader)

        row_list = list(tsv_reader)
        print('input length:',len(row_list))
        #print(row_list[0])

        for i in range(0,len(row_list)):
            row = row_list[i]
            converted_row = [int(ele) if ele.isdigit() else ele for ele in row] #convert int to int and str to str inside a row
            row_list[i] = converted_row
            #print(row)

        #print(row_list)
        
        if len(row_list) > 0:
            row = row_list[0]
            item = Item(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8])
            #print("item first:",item.feature)

        idx = 1
        while(idx <= len(row_list)):
            group = Group()
            group.circRNA = item

            #print(idx)
            if idx > len(row_list)-1 and item.feature == 'circRNA':
                exon_list = []
                group.exon_list = exon_list
                groups.append(group)
                break
            if idx > len(row_list)-1 and item.feature == 'exon':
                break

            row = row_list[idx]
            item = Item(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8])

            #print('dtype')
            #for val in row:
            #    print(type(val))

            exon_list = []
            cnt = idx+1
            #print("item:",item.feature)

            while item.feature == 'exon':
                exon_list.append(item)
                #print('idx=',idx)
                #print('cnt=',cnt)
                if cnt > len(row_list)-1:
                    break
                row = row_list[cnt]

                item = Item(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8])
                cnt = cnt + 1

            #print('exon list:',len(exon_list),'\n')
            #for i in range(0,len(exon_list)):
            #    exon = exon_list[i]
            #    print('exon:',exon.chr,exon.start,exon.end)
                
            group.exon_list = exon_list
            #print('\n')
            groups.append(group)
            idx = cnt

    return groups

#input = "simu_circ_star.gtf"
#output = "simu_circ_star_filtered.gtf"

input = sys.argv[1]
output = sys.argv[2]
min_coverage = float(sys.argv[3])

groups = input_to_groups(input)
print('len of groups in groups:',len(groups))

filtered_groups = []

for grp in groups:
    circRNA = grp.circRNA

    #attribute = circRNA.attribute
    #attribute_list = attribute.split(";")

    #if(len(attribute_list) < 3):
    #    continue
    
    #print(attribute_list[2])
    #cov_part = attribute_list[2]
    #cov_list = cov_part.split(" ")

    #if(len(cov_list) < 3):
    #    continue

    #quote_part = cov_list[2]
    #quote_list = quote_part.split("\"")

    #if(len(quote_list) < 3):
    #    continue
    
    #cov = int(quote_list[1])

    if(float(circRNA.score) > min_coverage):
        filtered_groups.append(grp)



print("size of filtered groups ",len(filtered_groups))
with open(output, "w") as f:
    for grp in filtered_groups:
        circRNA = grp.circRNA
        exon_list = grp.exon_list
        #print(len(exon_list))

        f.write(str(circRNA.chr) + "\t" + circRNA.src + "\t" + circRNA.feature + "\t" + str(circRNA.start) + "\t" + str(circRNA.end) + "\t" + str(circRNA.score) + "\t" + circRNA.strand + "\t" + circRNA.frame + "\t" + circRNA.attribute + "\n")

        for exon in exon_list:
            f.write(str(exon.chr) + "\t" + exon.src + "\t" + exon.feature + "\t" + str(exon.start) + "\t" + str(exon.end) + "\t" + str(exon.score) + "\t" + exon.strand + "\t" + exon.frame + "\t" + exon.attribute + "\n")

