#ifndef __potList_hh__
#define __potList_hh__

#include <cdsString.h>
#include <pot.h>

/**
 * Class which contains reference-counted Pots, 
 * with special methods for evaluating and reporting the potential terms.
 */
class PotList : public Pot {
public:
    PotList(const String name="");

    // adding/removing/listing potential terms in list
    void add( rc_Pot& pot);
    void remove(const String& instanceName);
    void removeAll();
    rc_Pot byName(const String& instanceName);
    CDSList< rc_Pot > list() { return list_; }

    //convenience (to make it act like a list)
    int size() const { return list_.size(); }
    rc_Pot operator[](int index) { return list_[index]; }
    const rc_Pot operator[](int index) const { return list_[index]; }

    // calculate energy, without and with derivs wrt atomic position
    // These also update the energyReports list
    //
    EnergyReport calcEnergy();
    EnergyReport calcEnergyAndDerivs(DerivList&);

    //methods to manipulate cached energy info
    //
    void addEnergyReport(const EnergyReport& v)    { energyReports_.append(v); }
    void clearEnergyReports()                      { energyReports_.resize(0); }
    CDSList<EnergyReport> energyReports() {return energyReports_;}

    String showReport() const;

private:
    CDSList< EnergyReport > energyReports_; 
    CDSList< rc_Pot > list_;
};

#endif /* __potList_hh__ */
