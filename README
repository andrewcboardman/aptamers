This program is intended to take sequencing data from an aptamer binding experiment (see Arter et al, Anal. Chem., 2018, 90 (17), pp 10302–10310), process it to identify sequences that successfully bound, and train an Ising model on these data. It consists of three parts:

Part I: Pre-processing
Trim adaptors, remove low quality sequences
One-hot encode
Standardise
Optionally: Augment FASTA with shifted sequences

Part II: Model training
This trains a model on the encoded data 
options: 
Independent sites
Mean-field ising
Convolutional neural net
Fully-connected neural net

Part III: Validation
This visualises the results of the model training
Generate samples from model (and compare to originals)
Plot inferred coefficients on histogram
Evaluate the probabilities of a set of sequences coming from the model


Prediction
This returns a suggested probability distribution to be used to generate new inputs (the marginal distribution for sequences that bind)


------------------------------------------------------------------------------

