##### ASSINGMENT #1 #####
# MAIN CODE

# Files required for each class
require './gene_Object.rb'
require './seed_stock_Object.rb'
require './hybrid_cross_Object.rb'

# If there are not 4 arguments in the command line a usage message is reported
if ARGV.length != 4
  abort("Usage: ruby process_database.rb  gene_information.tsv  seed_stock_data.tsv  cross_data.tsv  new_stock_file.tsv")
end

Gene.load_from_file(ARGV[0]) # Open file gene_information.tsv and use information of each gene as values of the class Gene properties
Seed_stock.load_from_file(ARGV[1]) # Open file seed_stock_data.tsv and use information of each gene in the seed stock as class Seed_stock attribute
Seed_stock.write_database(ARGV[3]) # Create file new_stock_file.tsv and writes updated information of seed stock after planting 
Hybrid_cross.load_from_file(ARGV[2]) # Open file cross_data.tsv and use information of each crossing as class Hybrid_cross attribute values
Hybrid_cross.chi_sq_test # Do a Chi-square test on the F2 cross data
linked=Gene.get_genes_linked # Return content of @linked_to attribute, as a hash with each gene id, if genetically-linked to another, and the id of the others it is linked to

puts;
puts "Final Report:"
puts;
# In case linked contains id values, print the names of genes linked. If not, print a message
if linked

  linked.each {|id_1, id_2| id_2.each { |id_2|puts "#{Gene.get_gene_info(id_1)[0]} is linked to #{Gene.get_gene_info(id_2)[0]}"}}
else
  puts "No gene is linked to another :("
end





