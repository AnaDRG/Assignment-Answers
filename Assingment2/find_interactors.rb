##### ASSINGMENT #2 #####
# Main script file

require './gene_Object.rb'
require './interaction_network_Object.rb'
require './annotation_Object.rb'
require 'net/http'

if ARGV.length != 2
  abort("Usage: ruby find_interactors.rb Genes_file.txt Report.txt")
end

def fetch(uri_str)  # "Fetch" routine that does some basic error-handling.
  begin
  address = URI(uri_str)  # create a "URI" object
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

def addtofile(file,protein,i) # Add to file information about each network interactor, where file = filename, protein=protein name and i=order in network
  gene=Gene.genes_in_file[Gene.prot_in_file.index(protein)] # Get Gene object from protein name
  file << "\tInteractor #{i}\n"
  file << "\tGene ID: #{gene.id}\n"
  file << "\tProtein Name: #{gene.prot}\n"
  file << "\tAnnotations:\n"
  gene.annot.each{|an|
    file << "\t\t#{an.type}: "
    file << "#{an.id}"
    (an.info.nil?) ? (file << "\n"):(file << "\t#{an.info}\n")
    }
  file << "\n" 
end

Gene.load_from_file(ARGV[0]) #ArabidopsisSubNetwork_GeneList.txt #ARGV[0] # Create Gene objects from genes in file
all_genes=Gene.genes_in_file # Get Gene objects just created
all_genes.each{|gene|  # For each Gene object, annotate them and find interactions with other proteins
  Gene.search_uniprot([gene,'g'])
  int=InteractionNetwork.new_network(gene.prot)
  }

new_file=ARGV[1] 
if File.exists?(new_file) # If the report file already exists, the program deletes it to update new_stock_file.tsv
      File.delete(new_file)
    end
new_f = File.new(new_file, "w") # Create file for writing
InteractionNetwork.in_file.each_with_index{|interaction,index| # interaction = Hash {root protein => [interactors]}
  proteins=interaction.values.join(", ")
  new_f << "Interaction #{index+1}\n#{interaction.keys[0]} interacts with #{proteins}\n"  
  interaction.each{|key,array_int| addtofile(new_f,key,1) and array_int.each_with_index do |interactor,i| # Add to file root protein and interactors in network information 
    addtofile(new_f,interactor,i+2)
  end
    }
}
new_f.close # Close file