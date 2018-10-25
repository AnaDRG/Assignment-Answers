##### ASSINGMENT #1 #####
# Gene class file

class Gene
  @@linked_number=0 # Set the number of genes that are linked to 0
  attr_accessor :id  # Gene id
  attr_accessor :name # Gene name 
  attr_accessor :mut_phen # Mutant phenotype
  attr_accessor :linked_to # Gene that is linked to
  
  def initialize (params = {}) # Get a value from the "new" call, or set a default
    # @id = "000000000", @name = "Unknown_name", @mut_phen = "Unknown_mutant_phenotype", @linked_to = "Not_linked_to_any_other_gene"
    @id = params.fetch(:id, "000000000")
    @name = params.fetch(:name, "Unknown_name") 
    @mut_phen = params.fetch(:mut_phen, "Unknown_mutant_phenotype")
    @linked_to = params.fetch(:linked_to, "Not_linked_to_any_other_gene")
  end
  
  def self.load_from_file(file) # Given a file name, it assigns each attribute value to file info
    
    abort("Sorry, I can't find #{file}") unless (File.exists?(file)) # Check if file exists 
    abort("Sorry, there's no info in #{file}") unless (File.size?(file)) # Check if file is empty (or presence, but that has already
    #been checked by empty?)
    f=File.open(file, "r") # Open file and read it
    f.readline # Header is removed from f    
    @id = Hash.new # Assings the property as a new hash
    @name = Array.new # Assings the property as a new array
    @mut_phen = Array.new
    
    f.each_line do |line| # For each line, assing values included to class properties
      line.chomp # Delete \n from lines
      next if line.strip.empty? # If an intermediate line is empty, pass to the next one
      id, name, mut_phen = line.split("\t") # Match each column of dataset
      id_format = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d/) # Check Arabidopsis gene identifier format and, if not right, returns a message and abort
      abort("Sorry, Arabidopsis gene identifiers have not the right format (e.g. AT1G12345)") unless id_format.match(id) 
      # Assing values in the id, name and mut_phen variables to class properties
      @id[id.to_s] = [name, mut_phen] # Id contains info about the other gene info 
      @name << name
      @mut_phen << mut_phen
      
    end
  end

  def self.get_gene_info(id) # Given a seed stock id, it returns all info about it from the gene data as an array
    return @id.fetch(id) if @id.key?(id)
  end
  
  def self.linked_to(id1,id2)
    
    @linked_to=Hash.new if @@linked_number==0 # If no linked genes have been identified, it assigns @linked_to as a new hash
    # If that gene is already linked to other genes, the new linked gene is added to its array associated @linked_to[id]
    # Genes linked are added to the id hash as an array: @id[ID][2] = array of ids from genes linked to the gene with the id ID
    (@linked_to[id1]) ? (@linked_to[id1] << [id2] and @id[id1][2].push([id2])):(@linked_to[id1]=[id2] and @id[id1].push([id2])) 
    (@linked_to[id2]) ? (@linked_to[id2] << [id1] and @id[id1][2].push([id2])):(@linked_to[id2]=[id1] and @id[id2].push([id1]))    
    @@linked_number+=1 # Count the amount of linked gene pairs
    
  end
  
  def self.get_genes_linked
    return @linked_to unless @@linked_number==0 # Get the hash linked_to only if there are genes linked (linked_number > 0)
  end
  
end
  
