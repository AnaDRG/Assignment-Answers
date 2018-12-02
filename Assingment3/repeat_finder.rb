 ## Assingment 3 ##
 
require 'net/http'
require 'bio'

def fetch(uri_str)  # "Fetch" routine that does some basic error-handling.
  begin
  address = URI(uri_str)  # Create a "URI" object
  response = Net::HTTP.get_response(address)  # Use the Net::HTTP object "get_response" method to call that address                                        
  case response   
  when Net::HTTPSuccess then  # When response Object is of type Net::HTTPSuccess
    return response  # Return that response object to the main code
  else
    raise Exception, "Something went wrong... the call to #{uri_str} failed; type #{response.class}"
    response = false
    return response  # Now we are returning 'False'
  end
  rescue
    abort("Sorry, something went wrong calling #{uri_str} :(")
  end
end

def load_from_file(file) # Given a file name, it assigns each attribute value to file info    
  abort("Sorry, I can't find #{file}") unless (File.exists?(file)) # Check if file exists 
  abort("Sorry, there's no info in #{file}") if (File.size?(file)==0) # Check if file is empty (or present, but that has already
  #been checked by empty?)
  all_genes=[] # Array of all genes in file
  f=File.open(file, "r") # Open file and read it    
  f.each_line do |id| 
    next if id.strip.empty? # If an intermediate line is empty, pass to the next one
    id_format = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d/) # Check Arabidopsis gene identifier format and, if not right, returns a message and abort
    abort("Sorry, Arabidopsis gene identifiers have not the right format (e.g. AT1G12345)") unless id_format.match(id)       
    all_genes.push(id.chomp)
  end
    return all_genes 
  end

def write_gff(file,parent,gene,c) # Write in file the location of the repetead motif, given the file name, the chromosome number,
  #the gene name and a bioseq feature that corresponds to the motif
  rep='Unknown_rep'
  strand='Unknown_strand'
  exon='Unkown_exon'
  c.qualifiers.each{|q| # Assign  the corresponding qualifier values of feature to rep (repeated motif) and strand (strand sign)
    rep=q.value if q.qualifier=='repeat_motif'
    strand=q.value if q.qualifier=='strand'
    exon=q.value if q.qualifier=='ex_id'
    }
    m=c.position.match(/(?<start>\d+)..(?<end>\d+)/) # Get coordinates of feature
    from=m[:start]
    to=m[:end]
    # Format GFF3: chromosome_number\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes
    file << "#{parent}\t.\trepeat_region\t#{from}\t#{to}\t.\t#{strand}\t.\tID=#{exon};parent=#{gene};target_sequence=#{rep}\n\n" 
end

def create_features(coord, bioseq, s, name, motif,exons) # Given an array of the repeated motif coordinates, a Bio::Sequence object bioseq and
    #the strand sign, name of repetition feature and motif, add a new feature to the sequence object 
    coord.each{|c| # For each coordinate, create a new feature where c is the location of motif 
    f = Bio::Feature.new(name, c ) 
    f.append(Bio::Feature::Qualifier.new('repeat_motif', motif))
    f.append(Bio::Feature::Qualifier.new('note', 'found by repeatfinder'))
    f.append(Bio::Feature::Qualifier.new('strand', s[coord.index(c)])) # Element in the strand sign array with index = coordinates c index
    f.append(Bio::Feature::Qualifier.new('ex_id', exons[coord.index(c)]))
    bioseq.features << f  # Append created feature
    } 
  
end

def create_gff(bioseq, file, rep, chr,gene) # Given a Bio::Sequence object bioseq, a file name, a feature name rep, the chromosome number
  #and gene name, call write_gff function for creating Summary file of the repeated motif
  bioseq.features.each do |feature|
      featuretype = feature.feature
      next unless featuretype == rep 
      write_gff(file,chr,gene,feature)
    end
end


if ARGV.empty? # At list, user must state gene list file
  abort("Usage: ruby repeat_finder.rb  gene_list.txt  [motif; default CTTCTT]")
end
gene_file = ARGV[0]
# If motif isn't set or isn't a nucleotide sequence, then set default cttctt 
(ARGV[1].nil? || (/[^atgc]/.match(ARGV[1].downcase) if ARGV[1])) ? motif = 'cttctt' : (motif = ARGV[1].downcase)
file_gff1= File.new("Summary_rep_gene_coord.gff", "w") # GFF3 file of motif coordinates in each gene
file_gff1 << "##gff-version 3\n"
file_gff2= File.new("Summary_rep_chr_coord.gff", "w") # GFF3 file of motif coordinates in chromosome
file_gff2 << "##gff-version 3\n"
file_nr=[]
ex_id=nil
exons=[]
all_genes=load_from_file(gene_file) # Return array of genes names from file
all_genes.each{|gene|  # For each gene, search for motif and if found, get coordinates and write them in file; if not found, then write gene
  #name in report file
  res = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}"); # Get EMBL file data of gene     
  record = res.body  # Get the "body" of the response     
  entry = Bio::EMBL.new(record) # Here we create an embl object  
  embl_bioseq = entry.to_biosequence # Convert a database entry to a Bio::Sequence object
  m=entry.accession.match(/chromosome:[A-Z0-9]+:(?<chromosome>\d+):(?<start_chr>\d+):(?<end_chr>\d+)/) #Get chromosome coordinates of gene
  chr= m[:chromosome]
  chr_start= m[:start_chr].to_i  
  chr_end= m[:end_chr].to_i
  coord_chr=[] # Array of motif coordinates in chromosome
  coord_gene=[] # Array of motif coordinates in gene
  exon_in_gene=[]
  sign=[]
  entry.features{|feature| # Search if feature corresponds to an exon
    if feature.feature='exon'
      unless /[A-Z]/.match(feature.position) # We are only interested in exons from gene (not remote entries)
        feature.qualifiers.each{|q| 
          if /exon_id=(.+)/.match(q.value)
            ex_id=$1
          else
            ex_id=nil
          end}
        next if ex_id.nil?
        begin
        exon=embl_bioseq.splicing(feature.position) #  If position of feature doesn't fit with gene sequence length, next feature 
        rescue
         next
        end
        exon_in_gene = feature.locations() # Bio::Locations object (container for Bio::Location objects)              
        exon.scan(/#{motif}/) do |repetition| # Search motif in the exon sequence
          from=$~.offset(0).first
          to=$~.offset(0).last
          
           if exon_in_gene.first.strand == -1 # Depending on strand sign, get coordinates
             length= to - from # Length of motif
             to=exon_in_gene.first.to - from + 1 # Where repetition ends in gene
             from= to - length + 1 # Where repetition starts in gene
             new_coord_chr= "#{from + chr_start - 2}..#{to + chr_start - 2 }" # Get coordinates in chromosome
             coord_chr << new_coord_chr
             new_coord_gene="#{from - 1}..#{to - 1}" # Get coordinates in gene
             coord_gene << new_coord_gene
             
           else
             new_coord_chr= "#{from + exon_in_gene.first.from + chr_start - 1}..#{to + chr_start + exon_in_gene.first.from - 2}" # Coordinates in gene 
             coord_chr << new_coord_chr
             new_coord_gene="#{from + exon_in_gene.first.from}..#{to + exon_in_gene.first.from - 1}" #Coordinates in chromosome
             coord_gene << new_coord_gene
             
           end
           s="++-"[exon_in_gene.first.strand <=> 0]
           sign << s
           exons << ex_id
        end
      end
    end
   }
  unless coord_chr.empty? # If motif not found, then add gene name to Report file; otherwise, write summary files
    #s="++-"[exon_in_gene.first.strand <=> 0] # Get strand sign where exon is located
    # Add new features to embl_bioseq corresponding to the repeated motif. Each feature corresponds to each motif found in sequence
    # Then create a summary file of the features just created
    # Follow the above steps for coordinates in gene
    create_features(coord_gene.uniq, embl_bioseq, sign, 'rep_in_gene', motif,exons) 
    create_gff(embl_bioseq,file_gff1, 'rep_in_gene', gene,gene)
    # Follow them for coordinates in chromosome
    create_features(coord_chr.uniq, embl_bioseq, sign, 'rep_in_chr', motif,exons)
    create_gff(embl_bioseq,file_gff2, 'rep_in_chr', chr,gene)
  else # Write a list of genes with no motif
    report="Report_no_motif.txt"
    unless File.exists?(report) # If the report file doesn't exists, create it  
      file_nr=File.new("Report_no_motif.txt", "w")
      file_nr<<"Genes without motif #{motif}:\n"      
    end
    file_nr << "#{gene}\n"    
  end
    
}
    
    
