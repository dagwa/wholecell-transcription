from roadrunner import RoadRunner, Logger

#Logger.setLevel(Logger.LOG_DEBUG)
Logger.setLevel(Logger.LOG_WARNING)

# Simulate model
rr = RoadRunner('/tmp/tx.sbml')
#rr = RoadRunner('/tmp/tx-reduced.sbml')
#rr = RoadRunner('/tmp/antsbml.xml')
results = rr.simulate(0, 3, 3, integrator='gillespie')

print(results)

with open('/tmp/sim_res.csv', 'w') as f:
    f.write(str(results))

# Plot
rr.plot()