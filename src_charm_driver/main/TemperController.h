#ifndef _TEMPERCONTROLLER_H_
#define _TEMPERCONTROLLER_H_
/*

   The Temper controller manages the temperature exchange process for tempers.

   Basic cycle is:
   if(iteration % TemperCycle==0)
   all tempers will pass their energies up to the TemperController
   all tempers will suspend
   - this is a global synchronization
   TemperController will engage the next temperature exchange.
   - Various schemes exist for this.
   - Configure choices wrt to temp exchange have no parallel driver impact.
   TemperController will issue new temperatures back to all Tempers.
   each temper will apply the new temperatures in Atoms and GSP
   each temper will resume
   -once all atoms and GSPs report completion in accepting the new temperature
   -Resumption within any temper is independent, so this is an
   instance level completion property.


 */

class TemperController : public CBase_TemperController
{
  public:
  TemperController(int simtype, double *temps, int tempsize, long seed, std::string history, std::string dirname);
    ~TemperController();
    TemperController(CkMigrateMessage *m) {}
    std::string history;
    std::string output_directory;
    void sumEnergies(EnergyStruct &inEnergy, int temper);
    void acceptData(int temper, EnergyStruct &energies);
    void acceptData();
    void output();
  private:
    int simType;
    int numTempers;
    int numBeads;
    int reportedIn;
    double *temperEnergy;
    double *temperatures;
    int *index;
    int switchdir;
    long seed;
    bool first;
};


#endif
