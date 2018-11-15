##### ASSINGMENT #2 #####
# Annotation class file

class Annotation
  attr_accessor :type # Type of info (i.e. GO_term)
  attr_accessor :id # Information ID (i.e. GO:0003006)
  attr_accessor :info # Information (i.e. developmental process involved in reproduction)  
  @@all_ids=[] # Class attribute, array with information IDs of all Annotation objects
  @@all_annotations=[] # Class attribute, array with all Annotation objects created
  
    def initialize (params = {}) # Get a value from the "new" call, or set a default
      # @id = "Unknown_gene_id", @id_prot = "Unknown_protein_id", @id_kegg = "Unknown_kegg_id"
      @type = params.fetch(:type, "Unknown_type")
      @id = params.fetch(:id, "Unknown_ID")
      @info = params.fetch(:info, "Unknown_information")
    end
    
    
    def self.annotate(*args) # Given a set of arguments, it returns an Annotation object
      t,i,inf=args # t= annotation type, i = annotation ID, inf = annotation information
      if Annotation.all_ids.include?(i) # If the ID is included in @@all_ids, then the annotation object already exists and is assigned to the variable an
        an=Annotation.all_annotations[Annotation.all_ids.index(i)]
      else # If the ID is not included in @@all_ids, then the annotation object is created and assigned to the variable an
        an=Annotation.new( :type => t, :id => i, :info => inf)
        @@all_ids << i # The ID of the new Annotation object is added to the @@all_id array
        @@all_annotations << an # The new Annotation object is added to the @@all_annotations array
      end
      return an
    end
    
    def self.all_ids # It returns the @@all_ids array
      return @@all_ids
    end
    def self.all_annotations # It returns the @@all_annotations array
      return @@all_annotations
    end
    
end

