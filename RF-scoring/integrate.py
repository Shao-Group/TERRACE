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


input1 = sys.argv[1] #our gtf file
input2 = sys.argv[2] #ML score file
output = sys.argv[3] #our gtf file with score embedded

groups1 = input_to_groups(input1)
print('len of groups in groups1:',len(groups1))

#print(groups1_dict.keys())

score_dict = {}
#read score file
with open(input2, "r", encoding="utf8") as score_file:
    tsv_reader = csv.reader(score_file, delimiter=",")

    # skip header
    next(tsv_reader)
    for row in tsv_reader:
        #print(row)
        (circRNA_id,score) = row
        #print(circRNA_id,score)
        score_dict.update({circRNA_id:float(score)})

with open(output, "w") as f:        
    for i in range(0,len(groups1)):
        group = groups1[i]
        circRNA = group.circRNA
        exon_list = group.exon_list

        hash = ""
        hash = hash + str(circRNA.chr) + ":" + str(int(circRNA.start)-1) + "|" + str(circRNA.end) + "|"
        for j in range(0,len(exon_list)):
            exon = exon_list[j]
            hash = hash + str(int(exon.start)-1) + "|" + str(exon.end) + "|"

        circRNA.score = score_dict[hash]

        f.write(str(circRNA.chr) + "\t" + circRNA.src + "\t" + circRNA.feature + "\t" + str(circRNA.start) + "\t" + str(circRNA.end) + "\t" + str(circRNA.score) + "\t" + circRNA.strand + "\t" + circRNA.frame + "\t" + circRNA.attribute + "\n")

        for exon in exon_list:

            exon.score = circRNA.score

            f.write(str(exon.chr) + "\t" + exon.src + "\t" + exon.feature + "\t" + str(exon.start) + "\t" + str(exon.end) + "\t" + str(exon.score) + "\t" + exon.strand + "\t" + exon.frame + "\t" + exon.attribute + "\n")


