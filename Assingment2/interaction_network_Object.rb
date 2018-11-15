##### ASSINGMENT #2 #####
# InteractionNetwork class file

class InteractionNetwork
  
  attr_accessor :members_l1  # Members of an interaction network at first level
  attr_accessor :members_l2 # Members of an interaction network at second level
  attr_accessor :root_protein # Members of an interaction network at third level
  @@in_file=[] # Class property, array with all protein interactors that are codified by genes in file
  
  def initialize (params = {}) # Get a value from the "new" call, or set a default
    # @membres_l1 = [], @members_l2 = [], @root_protein = "Unknown_root_protein"
    @members_l1 = params.fetch(:members_l1, [])
    @members_l2 = params.fetch(:members_l2, [])
    @root_protein = params.fetch(:root_protein, "Unknown_root_protein")
  end
 
  def self.load_interactors(prot_id)
    res = fetch("http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{prot_id}?format=tab25");  
    if res  # res is either the response object (Class Net::HTTPSuccess), or False, so you can test it with 'if' 
      body = res.body  # Get the "body" of the response
      members=Array.new
      # Only consider interactions between Arabidopsis thaliana proteins, with taxid:3702, repeated in text for each interaction four times in case both proteins are from A.thaliana)
      body.scan (/#{prot_id}|([A-Za-z0-9_]+)\(display_long\).+(psi-mi:"MI:[0-9]{4}").+((taxid:3702).+){4}.+intact-miscore:([0-9.]{4})/) do |i|
        unless i.nil?
          # i[0]= interactor name, i[1] = code of detection method, i.last = intact-miscore
          # Two hybrid array method is rejected because of the high number of FP and FN, with IDs MI:0397 (two hybrid array) and MI:0018 (two hybrid).
          # Score thereshold: 0.37 (typical value for two hybrid array)
          (i.last.to_f > 0.37 && i[1] != 'psi-mi:"MI:0397"'||i[1] !='psi-mi:"MI:0018"') ? (members << i[0] unless i[0].nil?) : ()
        end
      end
      (members.any?) ? (members=members.uniq and return members) : (return nil) # If not members return nil      
    else
      return nil
    end
  end
  
  def self.new_network(prot_id) # Get interactor members and assign as InteractionNetwork object properties  
      members = InteractionNetwork.load_interactors(prot_id) # Get level 1 members 
      included=Array.new # Array that includes proteins from the network in file
      if members
        int=InteractionNetwork.new( :root_protein => prot_id, :members_l1 => members) # If there are network members of level 1, create InteractionNetwork object
        members.map {|interactor|
          if Gene.prot_in_file.include?(interactor) # If interactor is in file, add to the included array
            included.push(interactor)            
          end         
          meml2=[]
          meml2 = InteractionNetwork.load_interactors(interactor) # Get level 2 members
          if meml2
            meml2.delete(prot_id) # Delete root protein from level 2 members
            int.members_l2.map { |array_ml2| # Delete members from l1 included in l2 interactions
            meml2.delete(int.members_l1[int.members_l2.index(array_ml2)]) if array_ml2.include?(interactor)
            }
          meml2.map{|mem| included.push(mem) if Gene.prot_in_file.include?(mem)} # If interactor is in file, add to the included array
          int.members_l2 << meml2         
          end
                   }
        unless included.empty? #If included not empty, add to @@in_file a hash with the root protein and the interactors from the network
          inc={int.root_protein => included.uniq}          
          @@in_file << inc
        end
        return int
      end     
  end 
  
  def self.in_file # Get the @@in_file array
    @@in_file.each{|i| i.each{|k,v| i.delete_if{|key,value|value == k and key == v}}} # Delete hash if repeated with a reverse order of members
    @@in_file.delete_if &:empty? # Delete empty hashes if there is any
    return @@in_file
  end
end


  
     














