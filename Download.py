import synapseclient 
import synapseutils
syn = synapseclient.Synapse() 
syn.login(authToken="")
files = synapseutils.syncFromSynapse(syn, "syn52663091", path="~/Desktop/NF") 