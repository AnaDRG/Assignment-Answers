
#### Assingment 4 ####
# Source of optimal e-value and coverage threshold:
# Moreno-Hagelsieb, G. and Latimer, K., 2007.
# Choosing BLAST options for better detection of orthologs as reciprocal best hits.
# Bioinformatics, 24(3), pp.319-324.

# Further analysis would require:
    # Search of GO terms: GO annotations for orthologous genes are likely to provide complementary information about conserved biological functions.
    # Build Clusters of Orthologous Groups (COGs) based on sequence similarity network of best reciprocal hits
    # Phylogenies and Tree Reconciliation.
    # Machine learning
    # Study conservation of local gene order (synteny), as a consequence of common ancestry that is most often observed among closely-related organisms.
    # Contrast results with available software and databases: Synergy, OrthoDB, Ensembl compara...


require 'bio'

abort("Usage: ruby orthologs_finder.rb  gene_database  protein_database evalue_exponent coverage") if ARGV[0] == "help"
(ARGV[0]) ? (file_nt=ARGV[0]) : (file_nt="TAIR10_seq_20110103_representative_gene_model_updated" and puts "Set TAIR10_seq_20110103_representative_gene_model_updated as gene database (default)")
(ARGV[1]) ? (file_aa=ARGV[1]) : (file_aa= "pep.fa" and puts "Set pep.fa as protein database (default)")
(ARGV[2] && ARGV[2].match(/[0-9.]+/) if ARGV[2]) ? (eval=ARGV[2]) : (eval= 6 and puts "Set evalue threshold as 1e-6 (default)")
(ARGV[3] && ARGV[3].match(/[0-9.]+/) if ARGV[3]) ? (cov=ARGV[3]) : (cov= 0.5 and puts "Set coverage threshold as 0.5 (default)")

abort("Sorry, I can't find #{file_nt}") unless (File.exists?(file_nt)) # Check if files exist
abort("Sorry, I can't find #{file_aa}") unless (File.exists?(file_aa)) 
abort("Sorry, there's no info in #{file_nt}") if (File.size?(file_nt)==0) # Check if files are empty (or present, but that has already
#been checked by empty?)
abort("Sorry, there's no info in #{file_aa}") if (File.size?(file_aa)==0)

system("makeblastdb -in #{file_nt} -dbtype 'nucl' -hash_index -out db_nt") # This creates some index files that are used by BLAST to speed-up the search
system("makeblastdb -in #{file_aa} -dbtype 'prot' -hash_index -out db_aa") 
fastacmd2 = Bio::Blast::Fastacmd.new("db_aa") # Retrieve FASTA formatted sequences from a blast database using NCBI fastacmd command
fastacmd1 = Bio::Blast::Fastacmd.new("db_nt")

factory2 = Bio::Blast.local('blastx', "db_aa", "-e #{eval}")   # The Bio::Blast class contains methods for running local or remote BLAST searches,
#as well as for parsing of the output of such BLASTs
factory1 = Bio::Blast.local('tblastn', "db_nt", "-e #{eval}") # blastx  compares the six-frame conceptual translation products of a nucleotide query sequence
#(both strands)against a protein sequence database while tblastn compares a protein query against the all six reading frames of a nucleotide sequence database.
$stderr.puts "Searching ... "  # Output to "standard error"

orthologs=Hash.new # Hash than contains the pairs of orthologs found
Bio::FlatFile.open(Bio::FastaFormat, file_aa) do |fasta_file| #  Automatically detect data format
  begin
   fasta_file.each do |entry| # For each protein in the protein file, do a blast search 
    begin
    report1 = factory1.query(">#{entry.entry_id}\n#{entry.seq}")
    hit1=report1.hits.first # Retrieve best hit
    # Get percent of the query sequence that overlaps the subject sequence: length of hit/ length of query sequence
    coverage1=(hit1.query_end.to_f + 1 - hit1.query_start.to_f)/report1.query_len.to_f 
    if coverage1 >= cov.to_f # Check if coverage under threshold
    fastacmd1.fetch(hit1.hit_id).each do |fasta| # Retrieve sequence of best hit from the database of gene sequences by hit ID
      begin
        report2=factory2.query(">\n#{fasta.seq}")
        hit2=report2.hits.first # Get best hit
         coverage2=(hit2.query_end.to_f + 1 - hit2.query_start.to_f)/report2.query_len.to_f # Calculate coverage
         if coverage2 >= cov.to_f and entry.definition == hit2.definition # Check if coverage under threshold and the orthology condition
          m=hit1.definition.match(/(?<gene>\S+)/) # Get gene ID
          gene=m[:gene]
          m=entry.definition.match(/(?<protein>[A-Za-z.0-9]+)/) # Get protein ID
          protein=m[:protein]
          orthologs[gene]=protein # Assing key as gene ID and value as protein ID if both orthologs
         end
      rescue
        next
      end
    end
    end
    rescue
      next
    end
   end
  rescue
    next
  end
end
unless orthologs.empty? # Write orthologs found in report
  report="Report_Orthologs.txt"
  file=File.new(report, "w")
  file << "Reciprocal Best Hits found (pairs of potential orthologs): #{orthologs.length}\n\nGENE\t\t\tPROTEIN\n"
  orthologs.each {|key,value| file << "#{key}\t\t#{value}\n"}       
  
end

