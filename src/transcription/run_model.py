from roadrunner import RoadRunner, Logger

Logger.setLevel(Logger.LOG_DEBUG)

# Simulate model
rr = RoadRunner('/tmp/tx.sbml')
results = rr.simulate(0, 10, integrator='gillespie')

# Plot
rr.plot()