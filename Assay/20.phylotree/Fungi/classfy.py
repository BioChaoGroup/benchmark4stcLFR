import sys
import os

levelDic = {0:"domain", 1:"superkingdom", 2:"kingdom", 3:"subkingdom", 4:"superphylum", 5:"phylum", 6:"subphylum", 7:"class", 8:"subclass", 9:"order", 10:"suborder", 11:"family", 12:"subfamily", 13:"genus", 14:"species", 15:"subspecies"}
#levelDic = {0:"domain", 1:"kingdom", 2:"phylum", 3:"class", 4:"order", 5:"family", 6:"genus", 7:"species", 8:"subspecies"}

levelDic_reverse = {}
for k, v in levelDic.items():
    levelDic_reverse[v] = k
    print(v, end="\t")
print("label")
#print(os.path.join(sys.path[0], 'log.txt'))
f = open(os.path.join(sys.path[0], 'log.txt'), 'w')

if __name__ == '__main__':
    test_file = r'paths.txt'
    with open(test_file, 'r') as fh:
        for line_index, line in enumerate(fh.readlines()):
            tmp = line.split(';')
            init_line = ['NA' for k in levelDic.keys()]
            last_index = 0
            for item in tmp[0:-1]:
                find = False
                item = item.strip()
                for key, values in levelDic_reverse.items():
                    if item.startswith(key):
                        init_line[values] = item
                        last_index = values+1
                        find = True
                        break
                if not find:
                    try:
                        init_line[last_index] = item
                        last_index += 1
                    except:
                        try:
                            f.write('Unknow class element in line: %s, colume: %s context: %s\n' % (line_index, last_index, item))
                        except:
                            pass

            init_line.append(tmp[-1].strip())
            print('\t'.join(init_line))
f.close()