import subprocess

def create_blast_db(fasta_file, db_name):
    subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl', '-out', db_name])

def run_blast(query, db, out_file):
    subprocess.run(['blastn', '-query', query, '-db', db, '-out', out_file, '-outfmt', '6 qseqid sseqid pident evalue'])

def parse_blast(file):
    best_hits = {}
    with open(file) as f:
        for line in f:
            columns = line.strip().split('\t')
            query_id = columns[0]
            subject_id = columns[1]
            identity = float(columns[2])
            evalue = float(columns[3])
            if query_id not in best_hits or evalue < best_hits[query_id]['evalue']:
                best_hits[query_id] = {'hit_id': subject_id, 'identity': identity, 'evalue': evalue}
    return best_hits

def find_reciprocal_best_hits(blast1, blast2):
    rbh = []
    for query_id, hit1 in blast1.items():
        hit_id = hit1['hit_id']
        if hit_id in blast2 and blast2[hit_id]['hit_id'] == query_id:
            rbh.append((query_id, hit_id, hit1['identity'], hit1['evalue'], blast2[hit_id]['identity'], blast2[hit_id]['evalue']))
    return rbh

def write_results(rbh, output_file):
    with open(output_file, 'w') as f:
        f.write("Query_ID\tHit_ID\tQuery_to_Hit_Identity\tQuery_to_Hit_E-value\tHit_to_Query_Identity\tHit_to_Query_E-value\n")
        for result in rbh:
            line_items = []
            for item in result:
                line_items.append(str(item))
            line = "\t".join(line_items) + "\n"
            f.write(line)


def main(file1, file2):
    db1 = "db1"
    db2 = "db2"
    out1 = "blast1.txt"
    out2 = "blast2.txt"

    create_blast_db(file1, db1)
    create_blast_db(file2, db2)

    print("Started Processing with Blast...")
    print("Please Wait...")
    run_blast(file1, db2, out1)
    print("Almost there :)...")
    run_blast(file2, db1, out2)

    blast1 = parse_blast(out1)
    blast2 = parse_blast(out2)

    rbh = find_reciprocal_best_hits(blast1, blast2)
    
    output_file = "reciprocal_best_hits.txt"
    write_results(rbh, output_file)
    print(f"RBH results written to {output_file}")

file1 = "E_coli_sRNAs.fasta"
file2 = "S_enterica_sRNAs.fasta"
main(file1, file2)