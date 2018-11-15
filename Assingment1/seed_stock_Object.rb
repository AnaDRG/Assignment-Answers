##### ASSINGMENT #1 #####
# Seed_stock class file

class Seed_stock

  attr_accessor :id  # Seed stock
  attr_accessor :mut_id # Mutant gene id - where we put Gene.id
  attr_accessor :l_planted # Last planted
  attr_accessor :storage # Storage name
  attr_accessor :grams # Grams remaining
  
  def initialize (params = {}) # Get a value from the "new" call, or set a default:
    #id = "000000000", mut_id = "Unknown_name", l_planted = "Unknown_mutant_phenotype", storage = "00/00/0000", grams = "0"
    @id = params.fetch(:id, "0000")
    @mut_id = params.fetch(:mut_id, "000000000") 
    @l_planted = params.fetch(:l_planted, "00/00/0000")
    @storage = params.fetch(:storage, "Unknown_storage")
    @grams = params.fetch(:grams, "0")
  end
  
  def self.load_from_file(file)
    abort("Sorry, I can't find #{file}") unless (File.exists?(file)) # Check if file exists 
    abort("Sorry, there's no info in #{file}") if (File.size?(file)==0) # Check if file is empty (or presence, but that has already
    #been checkd by empty?)
    f=File.open(file, "r") # Open file and read it
    $header = f.readline # Save header line    
    @id = Hash.new # Assings the property as a new hash
    @mut_id = Array.new # Assings the property as a new array
    @l_planted = Array.new
    @storage = Array.new
    @grams = Array.new    
    f.each_line do |line|
      
      line.chomp # Delete \n from lines
      next if line.strip.empty? # If an intermediate line is empty, pass to the next one
      id, mut_id, l_planted, storage, grams = line.split("\t") # Match each column of dataset
      @id[id.to_s] = [mut_id, l_planted, storage, grams.to_i] # Id contains info about the other info of seed stock
      # Indexes -> 0 = mut_id, 1 = l_planted, 2 = storage, 3 =grams
      @mut_id << mut_id
      @l_planted << l_planted
      @storage << storage
      @grams << grams.to_i
      
    end

  end

  def self.get_seed_stock(id) # Given a seed stock id, return all info from the data about it as an array
    return @id.fetch(id) if @id.key?(id)
  end
  
  
  def self.write_database (new_file)
    
    if File.exists?(new_file) # If the file already exists, the programme deletes it to update new_stock_file.tsv
      File.delete(new_file)
    end
    
    new_f = File.new(new_file, "w") # Create file for writing
    new_f << "#{$header}" # Write header (the same as the original file)
    t=Time.now.strftime("%d/%m/%Y") # Get actual date
    
    @id.each {|gene_id,stock_info| # For each id gene, update last time planted (stock_info[1]) and rest 7 grams of seeds (from stock_info[3]).
      # If there is not enough or 0, assing 0 gr of seeds and print message. 
      stock_info[1]=t
      (stock_info[3] - 7)<=0 ? (stock_info[3]=0 and puts "WARNING: we have run out of Seed Stock #{gene_id.to_s}"):(stock_info[3]-=7)
      # Print in file, for each row, the  gene id (gene_id) and the rest of info in stock_info (Mutant gene id, Last planted, Storage name 
      # and Grams remaining) separated by tab
      new_f << "#{gene_id}\t#{stock_info.join("\t")}\n"
      }
      
    new_f.close # Close file 
  end

end
