import sys
import os

#levelDic = {0:"domain", 1:"superkingdom", 2:"kingdom", 3:"subkingdom", 4:"superphylum", 5:"phylum", 6:"subphylum", 7:"class", 8:"subclass", 9:"order", 10:"suborder", 11:"family", 12:"subfamily", 13:"genus", 14:"species", 15:"subspecies"}
levelDic = {0:"domain", 1:"kingdom", 2:"phylum", 3:"class", 4:"order", 5:"family", 6:"genus", 7:"species", 8:"subspecies"}

levelDic_reverse = {}
for k, v in levelDic.items():
    levelDic_reverse[v] = k
    print(v, end="\t")
print("label")
#print(os.path.join(sys.path[0], 'log.txt'))
f = open(os.path.join(sys.path[0], 'log.txt'), 'w')


def filter_clade(value):
    if value.startswith('CLADE'):
        value = float(value.split('__').replace('CLADE1', '').replace('CLADE', ''))
        if value < 0.98:
            return True
    return False


def get_clade_level(value, value_map):
    if value.startswith('CLADE'):
        value = float(value.split('__')[0].replace('CLADE', ''))
        if value:
            for k, v in value_map.items():
                if v[0] <= value < v[1]:
                    return k, value
    return -1, 0


def get_value_map(level_dic):
    keys = list(levelDic.keys())
    sorted(keys)
    last_level = levelDic[keys[0]]
    last_key = keys[0]
    value_map = {}
    for k in keys[1:-1]:
        value_map[last_key] = [level_dic[last_level], level_dic[levelDic[k]]]
        last_level = levelDic[k]
        last_key = k
    value_map[last_key] = [level_dic[last_level], 1.1]
    return value_map


if __name__ == '__main__':
    test_file = r'paths.txt'
    DEFAULT_GE = 0.95
    DEFAULT_SP = 0.97
    test_file = sys.argv[1]
    if sys.argv[2] == "bac":
        DEFAULT_GE = 0.97
        DEFAULT_SP = 0.99
    else:
        DEFAULT_GE = 0.95
        DEFAULT_SP = 0.97
    claLevelDic = {"domain": 0.66, "kingdom": 0.7, "phylum": 0.76, "class": 0.77, "order": 0.83, "family": 0.89,
                   "genus": DEFAULT_GE, "species": DEFAULT_SP}
    value_map = get_value_map(claLevelDic)
    with open(test_file, 'r') as fh:
        for line_index, line in enumerate(fh.readlines()):
            tmp = line.split(';')
            init_line = ['NA' for k in levelDic.keys()]
            last_index = 0
            for item in tmp[0:-1]:
                find = False
                filtered = False
                item = item.strip()
                for key, values in levelDic_reverse.items():
                    # filtered = filter_clade(key)
                    # if filtered:
                    #     break
                    if item.startswith(key):
                        if init_line[values] != "NA":
                            f.write("Clade(%s) have in this place, line: %s, colume: %s context: %s, will replace, ignore\n" % (init_line[values], line_index, last_index, item))
                        init_line[values] = item
                        find = True
                        last_index = values+1
                        break
                if not find:
                    index, v = get_clade_level(item, value_map)
                    index1, v1 = get_clade_level(init_line[index], value_map)
                    if init_line[index] != "NA" and not init_line[index].startswith("CLADE"):
                        find = True
                    else:
                        if index != -1:
                            if index1 != -1 and v1 > v:
                                # origin value rate big than new value
                                pass
                            else:
                                init_line[index] = item
                                find = True
                                last_index = values + 1
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