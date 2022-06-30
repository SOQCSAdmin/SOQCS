import pysoqcs as soqcs

#
#   SOQCS TEMPLATE FILE
#
#


# Configure SOQCS
soqcs.cfg_soqcs(2)

# Create circuit and photon bunches
example = soqcs.qocircuit(2);
photons = soqcs.ph_bunch(example.num_levels(),1)
#Build a circuit

#
# BUILD YOUR CIRCUIT HERE
#    


# Create a simulator
simulator=soqcs.simulator()
# Simulate
measured=simulator.run(photons,example)

# Print measures
measured.show()