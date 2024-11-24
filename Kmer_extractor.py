import concurrent.futures
import os 
import time
import pandas
import concurrent
from Bio import SeqIO

def file_loader(filename,path):
    direction=path+filename
    sequence_cache=list(SeqIO.parse(direction,'fasta'))
    whole_genome=sequence_cache[0].seq
    organism=sequence_cache[0].id
    print(f'{filename} loaded')
    return whole_genome,organism

def kmer_extractor(sequence,kmer,size,organism):
    def kmer_segmentation(sequence,kmer,size):
        print('{0:-^50}'.format('Kmer extraction'))
        print('extracting kmers')
        kmerized_genome=[]
        for i in range(size):
            fragment=sequence[i:i+kmer]
            kmerized_genome.append(fragment)
        print('kmers extracted')
        print('{0:-^50}'.format(''))
        print('')
        print('{0:-^25}'.format(''))
        print('sorting kmers')
        kmerized_genome.sort()
        print('kmers sorted')
        print('{0:-^25}'.format(''))
        print('')
        return kmerized_genome

    def kmer_processing(kmerized_genome):
        print('{0:-^50}'.format('Searchlist creation'))
        print('creating search list')
        processing_frame=pandas.DataFrame({'col':kmerized_genome})
        processing_frame.drop_duplicates(inplace=True)
        search_list=processing_frame['col'].tolist()
        print('search list created')
        print('{0:-^50}'.format(''))
        print('')
        return search_list
    
    def kmer_search(kmerized_genome,target,size,organism):
        low=0
        high=size

        def First_occurrence(kmerized_genome,low,high,target):
            if high>=low:
                mid=(low+high)//2

                if (mid==0 or target>kmerized_genome[mid-1]) and kmerized_genome[mid]==target:
                    return mid
                elif target>kmerized_genome[mid]:
                    return First_occurrence(kmerized_genome,mid+1,high,target)
                else:
                    return First_occurrence(kmerized_genome,low,mid-1,target)
            return -1

        def Last_occurrence(kmerized_genome,size,low,high,target):
            if high>=low:
                mid=(low+high)//2

                if (mid==size-1 or target<kmerized_genome[mid+1]) and kmerized_genome[mid]==target:
                    return mid
                elif target<kmerized_genome[mid]:
                    return Last_occurrence(kmerized_genome,size,low,mid-1,target)
                else:
                    return Last_occurrence(kmerized_genome,size,mid+1,high,target)
            return -1

        first=First_occurrence(kmerized_genome,low,high,target)
        last=Last_occurrence(kmerized_genome,size,low,high,target)
        result=(last-first)+1
        return result

    segmented_genome=kmer_segmentation(sequence,kmer,size)
    processed_genome=kmer_processing(segmented_genome)
    count_cache=[]
    filename='K-merized-'+organism+'-'+'.csv'
    print('{0:-^50}'.format('Target counting'))
    print('searching targets')
    for i in range(len(processed_genome)):
        result=kmer_search(segmented_genome,processed_genome[i],size,organism)
        count_cache.append(result)
    print('targets found')
    print('{0:-^50}'.format(''))
    print('')
    print('{0:-^50}'.format('Output file'))
    print('saving output file')
    report=pandas.DataFrame({'kmers':processed_genome,organism:count_cache})
    report.to_csv(filename)
    print('output file saved')
    print('{0:-^50}'.format(''))
    
def main():
    os.system('clear')
    begin=time.perf_counter()
    path='/home/renegade/Documents/MCBCI/Tesis-Axel/Software/Kmer_extractor/bacterial_genomes/' #cambiar al path que vayan a usar 
    available_files=os.listdir(path)
    kmer=13
    path_list=[]
    kmer_list=[]
    n_list=[]
    genomes=[]
    organisms=[]
    
    for _ in range(len(available_files)):
        path_list.append(path)
        kmer_list.append(kmer)

    print('{0:-^50}'.format('Loaded files'))
    with concurrent.futures.ProcessPoolExecutor() as overlord:
        results=list(overlord.map(file_loader,available_files,path_list))
    print('{0:-^50}'.format(''))

    for _ in range(len(results)):
        temp_index=results[_]
        n_list.append(len(temp_index[0]))
        genomes.append(temp_index[0])
        organisms.append(temp_index[1])
    
    with concurrent.futures.ProcessPoolExecutor() as overlord:
       overlord.map(kmer_extractor,genomes,kmer_list,n_list,organisms)

    end=time.perf_counter()
    print("")
    print(f"Program finished in {end-begin} seconds")

main()
