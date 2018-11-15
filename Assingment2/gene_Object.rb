##### ASSINGMENT #2 #####
# Gene class file

class Gene
  
  attr_accessor :id  # Gene id
  attr_accessor :prot # Protein name
  attr_accessor :annot # Annotation objects associated to gene
  attr_accessor :kegg_entry # Kegg entry ID, unique for each Gene object  
  @@all_proteins=[] # Class attribute, array with all proteins codified by file genes
  @@all_genes=[] # Class attribute, array with all gene objects from file genes
  
  def initialize (params = {}) # Get a value from the "new" call, or set a default
    # @id = "Unknown_gene_id", @prot="Unknown_protein_name", @annot="Unknown_annotation", @kegg_entryº=nil
    @id = params.fetch(:id, "Unknown_gene_id")
    @prot = params.fetch(:prot, "Unknown_protein_name")
    @annot = params.fetch(:annot, "Unknown_annotation")
    @kegg_entry = params.fetch(:kegg_entry, nil)
  end
  
  def self.load_from_file(file) # Given a file name, it assigns each attribute value to file info    
    abort("Sorry, I can't find #{file}") unless (File.exists?(file)) # Check if file exists 
    abort("Sorry, there's no info in #{file}") if (File.size?(file)==0) # Check if file is empty (or present, but that has already
    #been checked by empty?)
    #all_genes=[]
    f=File.open(file, "r") # Open file and read it    
    f.each_line do |id| # For each line, assing values included to class properties
      next if id.strip.empty? # If an intermediate line is empty, pass to the next one
      id_format = Regexp.new(/A[Tt]\d[Gg]\d\d\d\d\d/) # Check Arabidopsis gene identifier format and, if not right, returns a message and abort
      abort("Sorry, Arabidopsis gene identifiers have not the right format (e.g. AT1G12345)") unless id_format.match(id)
      # Create a new gene object and assign ID without \n with .chomp 
      g=Gene.new(:id => id.chomp)     
      @@all_genes.push(g)
    end    
  end
    
  def self.search_kegg(k_entry) # Search KEGG Pathway and return pathway ID and pathway name
      res = fetch("http://togows.org/entry/kegg-compound/#{k_entry}") unless k_entry.nil?;
      if res  # res is either the response object (Class Net::HTTPSuccess), or False, so you can test it with 'if'
        body = res.body  # get the "body" of the response
        m=body.match(/PATHWAY\s+(?<path_id>[a-zA-Z0-9]+)\s+(?<path_name>[a-zA-Z0-9 -]+)\n/)                
        (m) ? (return m[:path_id], m[:path_name]):(return "Unknown_pathway") # Return pathway ID and pathway name or "Unknown_pathway"    
      else
        return "Unknown_pathway"
      end
    end
  
  def self.search_uniprot(gp) # Search protein name knowing gene ID (or vice versa), kegg entry ID and GO terms (biological process). gp= [gene_ID/protein_name, g (gene)/p (protein)]
      id=gp[0].id if gp[1] == 'g'
      id=gp[0] if gp[1] == 'p'
      res = fetch("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=uniprotkb&id=#{id}&format=default&style=raw");
      if res  # res is either the response object (Class Net::HTTPSuccess), or False, so you can test it with 'if'
        body = res.body  # get the "body" of the response
        if gp[1] =='g' # Search protein name and add as an object property
          m=body.match(/ID\s+(?<prot_name>[a-zA-Z0-9_]+)/)
          (m) ? (gp[0].prot = m[:prot_name].downcase and @@all_proteins << m[:prot_name].downcase) :(@@all_proteins << "Unknown_protein")
          g=gp[0]
        else # Search gene ID, create Gene object and add id and prot object properties
          m=body.match(/OrderedLocusNames=(?<gene_id>[a-zA-Z0-9_]+)/)
          g=Gene.new(:id => m[:gene_id], :prot => id.downcase) if m
        end
        # Look for KEGG entry code needed for searching interactions
        m=body.match(/KEGG;\s(?<k_entry>[a-zA-Z0-9_:]+)/)
        if m                    
          g.kegg_entry=m[:k_entry] # Assign KEGG entry as Gene object property
          id_path, name_path=Gene.search_kegg(m[:k_entry]) # Get pathway information from KEGG entry
          # Add to annot property an Annotation Object created from pathway information. If array is not created, create it with this object.
          g.annot << Annotation.annotate("Pathway",id_path, name_path) rescue g.annot=[Annotation.annotate("Pathway",id_path, name_path)]        
        end
        # Search GO terms from the biological_process part of the GO Ontology
        go=body.scan(/GO;\s(GO:[0-9]{7});\sP:([A-Za-z0-9 -]+);/) # go = array where each element is a hash {GO ID => GO term}
        if go.any? # If there are GO terms associated, add to annot property an Annotation Object created from GO term information. If array is not created, create it with this object.
          go=go.uniq
          go.each {|go| g.annot << Annotation.annotate("GO_term_Biol_process", go[0], go[1]) rescue g.annot=[Annotation.annotate("GO_term_Biol_process",go[0], go[1])]} 
        end
        return g if gp[1] == 'p' 
      end      
  end
    
  def self.annotate_interactors(id) # Use protein names from interactors (not in file) to get information about them (not used in main script but useful for annotation of the full networks obtained)
      Gene.search_uniprot([id,'p'])
  end

  def self.genes_in_file # Get @@all_genes, array with all Gene object created from file  
      return @@all_genes
    end

  def self.prot_in_file # Get @@all_genes, array with all Gene object created from file  
      return @@all_proteins
    end
  
end