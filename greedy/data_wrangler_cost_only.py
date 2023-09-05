import sys
import subprocess


def construct_graph(backbone_seq, start_pos, end_pos, pos, ref, alts, graph_file_name):

    # add backbone edges
    graph_file = open(graph_file_name, 'w')
    for i in range(len(backbone_seq)):
        graph_file.write(str(int(start_pos)+i) + ' ' + str(int(start_pos)+i+1) + ' ' + backbone_seq[i] + ' ' + '-' + '\n')

    # add alt paths
    new_vertex = int(end_pos) + 2
    for i in range(len(pos)):

        for a in alts[i]:
            for j, sym in enumerate(a):
                if j == 0:
                    start = pos[i]

                else:
                    start = str(new_vertex)

                if j == len(a)-1:
                    end = str(int(pos[i]) + len(ref[i]))
                else:
                    new_vertex += 1
                    end = str(new_vertex)

                graph_file.write(start + ' ' + end + ' ' + sym + ' ' + str(i) + '\n')

    graph_file.close()


def extract_variants(vcf_file_name, chr, start_pos, end_pos):

    # bcftools='/storage/coda1/p-saluru8/0/ntavakoli6/hged/software/bcftools-1.9/bcftools'
    bcftools='/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools'
    p1 = subprocess.Popen([bcftools, 'query',
                           '--include', '(TYPE="snp" || TYPE="indel") && GT="alt"',
                           '--regions', chr+':'+start_pos+'-'+end_pos,
                           '--format', '[%POS %REF %ALT %SAMPLE %GT *]', vcf_file_name], stdout=subprocess.PIPE)
    stdout = p1.communicate()
    lines = stdout[0].decode("utf-8").split('*')
    pos = ['']*(len(lines)-1)
    ref = ['']*(len(lines)-1)
    alt = ['']*(len(lines)-1)
    sample = ['']*(len(lines)-1)
    haplotype = ['']*(len(lines)-1)
    for i, line in enumerate(lines):
        line.strip()
        if len(line) > 0:
            pos[i], ref[i], alt[i], sample[i], haplotype[i] = line.strip().split(' ')

    return pos, ref, alt, sample, haplotype


def obtain_substrings(backbone_seq, start_pos, alpha, pos, sample_alts, pos_substring_file_name):

    f = open(pos_substring_file_name, 'w')
    final_pos = int(start_pos) + len(backbone_seq)
    start_pos = int(start_pos)
    print()
    for i, p in enumerate(pos):

        string_set = {}
        for tup in sample_alts[p]:
            sample = tup[0]
            for gt in {0, 1}:
                if (gt == 0 and tup[1][0] != '0') or (gt == 1 and tup[1][2] != '0'):

                    # obtain substring
                    num_variant_added = 0
                    s = ''
                    current_pos = int(p)
                    while current_pos < final_pos and len(s) < alpha:
                        if str(current_pos) in sample_alts:
                            sample_found_at_pos = False
                            for tup2 in sample_alts[str(current_pos)]:
                                if tup2[0] == sample:
                                    if gt == 0:
                                        alt_choice = int(tup2[1][0])
                                    if gt == 1:
                                        alt_choice = int(tup2[1][2])
                                    if alt_choice != 0:
                                        s += tup2[3].split(',')[alt_choice-1]
                                        current_pos += len(tup2[2]) # add reference length
                                        sample_found_at_pos = True
                                        num_variant_added += 1

                            if not sample_found_at_pos:
                                s += backbone_seq[current_pos - start_pos]
                                current_pos += 1
                        else:
                            s += backbone_seq[current_pos - start_pos]
                            current_pos += 1


                    # add to this positions set
                    if s not in string_set:
                        string_set[s] = 1

        f.write(p + ' ' + ' '.join(string_set) + '\n')
    f.close()


def construct_pos_ref_alt(pos, ref, alt):
    n = len(pos)
    pos_new = [pos[0]]
    ref_new = [ref[0]]
    alt_new = [alt[0].split(',')]
    for i in range(1, n):
        if pos[i] == pos[i-1] and ref[i] == ref[i-1] and alt[i] == alt[i-1]:
            continue
        else:
            pos_new.append(pos[i])
            ref_new.append(ref[i])
            alt_new.append(alt[i].split(','))
    return pos_new, ref_new, alt_new


def construct_pos_sample_alts(pos_raw, ref_raw, alts_raw, samples_raw, gt_raw, max_num_variants):

    sample_alt = {}
    for i in range(len(pos_raw)):
        if pos_raw[i] in sample_alt:
            sample_alt[pos_raw[i]].append((samples_raw[i], gt_raw[i], ref_raw[i], alts_raw[i]))
        else:
            sample_alt[pos_raw[i]] = [(samples_raw[i], gt_raw[i], ref_raw[i], alts_raw[i])]
        if len(sample_alt) == max_num_variants:
            return sample_alt
    return sample_alt


def get_backbone(reference_file_name, chr, start_pos, end_pos):

    # samtools='/storage/coda1/p-saluru8/0/ntavakoli6/hged/software//samtools-1.12/samtools'
    # bcftools='/storage/coda1/p-saluru8/0/ntavakoli6/hged/software/bcftools-1.9/bcftools'
    samtools='/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools'
    bcftools='/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools'

    p1 = subprocess.Popen([samtools, 'faidx', reference_file_name, chr + ':' + start_pos + '-' + end_pos],
                          stdout=subprocess.PIPE)
    stdout = p1.communicate()
    lines = stdout[0].decode("utf-8").split('\n')
    backbone_seq = ''
    for line in lines[1:]:
        line = line.strip()
        backbone_seq += line

    return backbone_seq


def get_cost_values(pos, ref, alts, greedy_cost_file):

    f = open(greedy_cost_file ,'w')
    for i, p in enumerate(pos):
        delta = max([max(len(a), len(ref[i])) for a in alts[i]])
        f.write(p + ' ' + str(delta) + '\n')
    f.close()


if __name__ == "__main__":

    reference_file_name = sys.argv[1]
    vcf_file_name = sys.argv[2]
    chr = sys.argv[3]
    start_pos = sys.argv[4]
    end_pos = sys.argv[5]
    max_num_variants = int(sys.argv[6])
    alpha = int(sys.argv[7])
    graph_out_file_name = sys.argv[8]
    # pos_substring_file_name = sys.argv[9]
    greedy_cost_file_name = sys.argv[9]

    print('Initial query')
    pos_raw, ref_raw, alts_raw, samples_raw, gt_raw = extract_variants(vcf_file_name, chr, start_pos, end_pos)

    print('Getting pos, ref, alt')
    pos, ref, alts = construct_pos_ref_alt(pos_raw, ref_raw, alts_raw)
    pos = pos[:max_num_variants]
    print('Number positions: ', len(pos))

    # print('Getting sample table')
    # sample_alts = construct_pos_sample_alts(pos_raw, ref_raw, alts_raw, samples_raw, gt_raw, max_num_variants)

    # # get backbone string
    print('Getting backbone')
    end_pos = str(int(end_pos) + len(ref[-1]))
    backbone_seq = get_backbone(reference_file_name, chr, start_pos, end_pos)

    # # write file with pos_substrings
    # print('Obtaining substrings')
    # # obtain_substrings(backbone_seq, start_pos, alpha, pos, sample_alts, pos_substring_file_name)

    # construct graph file
    print('Obtaining graph file')
    construct_graph(backbone_seq, start_pos, end_pos, pos, ref, alts, graph_out_file_name)

    # get delta's for greedy algorithm
    print('Writing cost file for greedy algorithm')
    get_cost_values(pos, ref, alts, greedy_cost_file_name)
