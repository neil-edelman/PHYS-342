struct Simulation;

struct Simulation *Simulation(int size, void (*v)(float, float, float));
void Simulation_(struct Simulation **simulationPtr);
int SimulationGetSize(const struct Simulation *simulation);
float SimulationGetMu(const struct Simulation *s, const int x, const int y, const int z);
int (*SimulationGetExplode(const struct Simulation *s))(struct Simulation *, const int);
void SimulationClearExplode(struct Simulation *s);
int SimulationUpdate(struct Simulation *s);
int SimulationCurrent(const struct Simulation *s);
int SimulationMagnetic(const struct Simulation *s);
void SimulationAnimation(struct Simulation *s);
