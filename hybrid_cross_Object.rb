##### ASSINGMENT #1 #####
# Hybrid_cross class file

class Hybrid_cross

  attr_accessor :p1p2  # Parent1 and Parent2 seed stock ids
  attr_accessor :f2_w # Number of F2 wildtype obtained
  attr_accessor :f2_p1 # Number of F2 with parent 1 genotype obtained
  attr_accessor :f2_p2 # Number of F2 with parent 2 genotype obtained
  attr_accessor :f2_p1p2 # Number of F2 with both parents genotype obtained
  attr_accessor :linked # Number of F2 with both parents genotype obtained
  
  def initialize (params = {}) # Get a value from the "new" call, or set a default:
    
    #p1p2 = ["Unknown parent 1", "Unknown parent 2"], f2_w = 0, f2_p1 = 0, f2_p2 = 0, f2_p1p2 = 0, linked = "No genes linked"  
    @p1p2 = params.fetch(:p1p2, ["Unknown parent 1", "Unknown parent 2"] )
    @f2_w = params.fetch(:f2_w, 0)
    @f2_p1 = params.fetch(:f2_p1, 0)
    @f2_p2 = params.fetch(:f2_p2, 0)
    @f2_p1p2 = params.fetch(:f2_p1p2, 0)
    @linked = params.fetch(:f2_p1p2, "No genes linked")
    
  end
  
  def self.load_from_file(file)
    
    abort("Sorry, I can't find #{file}") unless (File.exists?(file)) # Check if file exists 
    abort("Sorry, there's no info in #{file}") unless (File.size?(file)) # Check if file is empty (or presence, but that has already
    #been checkd by empty?)
    f=File.open(file, "r") # Open file and read it
    f.readline # Header is removed from f    
    @p1p2 = Hash.new # Assings the property as a new hash
    @f2_w = Array.new # Assings the property as a new array
    @f2_p1 = Array.new
    @f2_p2 = Array.new
    @f2_p1p2 = Array.new
    @linked = Hash.new    
    f.each_line do |line|
      
      line.chomp # Delete \n from lines
      next if line.strip.empty? # If an intermediate line is empty, pass to the next one
      p1, p2, f2_w, f2_p1, f2_p2, f2_p1p2 = line.split("\t") #match each column of dataset      
      @p1p2[[p1.to_s,p2.to_s]] = [f2_w.to_i, f2_p1.to_i, f2_p2.to_i, f2_p1p2.to_i] # Ids of the parents allow to know the number of each
      # F2 genotype. 
      @f2_w << f2_w.to_i
      @f2_p1 << f2_p1.to_i
      @f2_p2 << f2_p2.to_i
      @f2_p1p2 << f2_p1p2.to_i
      
    end    
  end
    
    def self.chi_sq_test()
      
      freq=['9/16'.to_r, '3/16'.to_r, '3/16'.to_r, '1/16'.to_r] # Frequences expected if genes not linked      
      @p1p2.each {|parents_id, genotypes| # Key = each pair of gene ids of parent 1 and parent 2 in an array [id_parent1, id_parent2]; Value =
        # amount of all kind of F2 genotypes in an array [f2_w, f2_p1, f2_p2, f2_p1p2]
        chisq=0 # For each pair of gene "parents", assing an initial value 0
        genotypes.each_index {|i| # For each genotype, compute the chi square values, as the summatory of (observed - expected)^2 / expected
           # where observed = the amount of each genotype (ie f2_w) and expected = the expected frequency of that genotype (in the array freq)
           # multiplied by the total amount of F2 genotypes observed (the sum of all values in the array of genotypes = each "genotypes" array per
           # pair of gene parents = each "parents_id" array)
          chisq += ((genotypes[i] - (freq[i] * genotypes.inject(:+)))**2 / (freq[i] * genotypes.inject(:+)))           
          }
        if (chisq >= 8) then # Value of chi square needed to refuse null hypothesis (genes not linked) for 3 grades of freedom
          
          begin # In case of error, becasuse some info needed is missed, a message appears and abort (rescue block)
          g_id1=Seed_stock.get_seed_stock(parents_id[0]) # Knowing the gene id (parents_id[index]), get the info about that gene (mutant id included and saved in g_id)
          g_id2=Seed_stock.get_seed_stock(parents_id[1]) # Parent 1: parents_id [0] and g_id1; Parent 2: parents_id [1] and g_id2
          p1=Gene.get_gene_info(g_id1[0])[0] # Knowing the mutant id (g_id[0]), get the info about that gene (gene name included and saved in p1 and p2)  
          p2=Gene.get_gene_info(g_id2[0])[0] #Parent 1 gene name: p1; Parent 2 gene name: p2          
          puts "Recording: #{p1} is genetically linked to #{p2} with chisquare score #{chisq.to_f}"
          @linked[[p1,p2]] = chisq.to_f # Assinged to linked property for the names of genes linked the value of chi-square
          Gene.linked_to(g_id1[0],g_id2[0]) # Relate both gene linked by ids (mutant ids)
          rescue
            abort("There is not enough data for full genetic linkage analysis, sorry :( Check there is info of all genes in the three files used")
          end     
            
        end        
        }
    end

  end
