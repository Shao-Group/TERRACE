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

#input1 = "simu_sorted.gtf"
#input2 = "simu_circ_star_sorted.gtf"

input1 = sys.argv[1]
input2 = sys.argv[2]
#output = sys.argv[3]

groups1 = input_to_groups(input1)
print('len of groups in groups1:',len(groups1))

#print('group1 info.........')
#for i in range(0,len(groups1)):
#    group = groups1[i]
#    print('circRNA:',group.circRNA.chr,group.circRNA.start,group.circRNA.end)
#    print('exons:')
#    print(len(group.exon_list))
#
#    for j in range(0,len(group.exon_list)):
#        exon = group.exon_list[j]
#        print('exon:',exon.chr,exon.start,exon.end)

groups2 = input_to_groups(input2)
print('len of groups in groups2:',len(groups2))

# print('group2 info.........')
# for i in range(0,len(groups2)):
#     group = groups2[i]
#     print('circRNA:',group.circRNA.chr,group.circRNA.start,group.circRNA.end)
#     print('exons:')
#     print(len(group.exon_list))

#     for j in range(0,len(group.exon_list)):
#         exon = group.exon_list[j]
#         print('exon:',exon.chr,exon.start,exon.end)

groups1_dict = {}

for i in range(0,len(groups1)):
    group = groups1[i]
    circRNA = group.circRNA
    exon_list = group.exon_list

    hash = ""
    hash = hash + str(circRNA.chr) + "|" + str(circRNA.start) + "|" + str(circRNA.end) + "|"
    for j in range(0,len(exon_list)):
        exon = exon_list[j]
        hash = hash + str(exon.chr) + "|" + str(exon.start) + "|" + str(exon.end) + "|"

    if hash in groups1_dict.keys():
        lst = list(groups1_dict[hash])
        lst[1] = lst[1] + 1
        groups1_dict[hash] = tuple(lst)
    else:
        groups1_dict.update({hash:(group,1)})
    #print (hash)

#print(groups1_dict)

groups2_dict = {}

for i in range(0,len(groups2)):
    group = groups2[i]
    circRNA = group.circRNA
    exon_list = group.exon_list
    #print(len(exon_list))

    hash = ""
    hash = hash + str(circRNA.chr) + "|" + str(circRNA.start) + "|" + str(circRNA.end) + "|"
    for j in range(0,len(exon_list)):
        exon = exon_list[j]
        hash = hash + str(exon.chr) + "|" + str(exon.start) + "|" + str(exon.end) + "|"

    if hash in groups2_dict.keys():
        lst = list(groups2_dict[hash])
        lst[1] = lst[1] + 1
        groups2_dict[hash] = tuple(lst)
    else:
        groups2_dict.update({hash:(group,1)})
    #print (hash)

#print(groups2_dict)

# for ky in groups2_dict.keys():
#     prinet(len(groups2_dict[key][0].exon_list))

print("length of dict/distinct circRNA groups1:",len(groups1_dict))
print("length of dict/distinct circRNA groups2:",len(groups2_dict))

#find similarity and diff between two dicts

dict1_set = set(groups1_dict)
dict2_set = set(groups2_dict)

match = dict1_set.intersection(dict2_set)
#match = groups1_dict.keys() & groups2_dict.keys()

match_len = len(match)
print("matching items:",match_len)

match_dict = dict.fromkeys(match, 1)
#print(len(match_dict))

# for val in match:
#     print(val,groups1_dict[val])

circRNA_count_groups1 = len(groups1_dict)
circRNA_count_groups2 = len(groups2_dict)

match_gr1 = match_len/circRNA_count_groups1*100
print("Percentage of match against group 1/Recall:",match_gr1)

match_gr2 = match_len/circRNA_count_groups2*100
print("Percentage of match against group 2/Precision:",match_gr2)