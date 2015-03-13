from roadrunner import RoadRunner, Logger

#Logger.setLevel(Logger.LOG_DEBUG) # too verbose
Logger.setLevel(Logger.LOG_WARNING)

# Load model

rr = RoadRunner('/tmp/tx.sbml')
#rr = RoadRunner('/tmp/tx-reduced.sbml')
#rr = RoadRunner('/tmp/antsbml.xml')

# Simulate model
results = rr.simulate(0, 10, 30, integrator='gillespie')

print(results)

# Write results to temp file
with open('/tmp/sim_res.csv', 'w') as f:
    f.write(str(results))

# Plot
rr.plot()