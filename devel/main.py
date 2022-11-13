import pysoqcs as soqcs

#
#   SOQCS TEMPLATE FILE
#
#

# Create circuit
example = soqcs.qodev(2,2);


#Build a circuit
#
# BUILD YOUR CIRCUIT HERE
#    


# Create a simulator
simulator=soqcs.simulator()
# Simulate
measured=simulator.run(example)

# Print measures
measured.show()