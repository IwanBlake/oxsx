#include <catch.hpp>
#include <ROOTNtuple.h>
#include <IO.h>
#include <OXSXDataSet.h>
#include <iostream>
TEST_CASE("Writing a data set to disk  and reading back"){
    
    // Create fake data set
    OXSXDataSet origDataSet;
    std::vector<double> eventObs(4);

    for(unsigned i = 0; i < 2181126; i++){
        if(!(i%1000000))
            std::cout << 100. * i/2181126 << " % " << std::endl;
        for(size_t j = 0; j < eventObs.size(); j++)
            eventObs[j] = j;
        
        origDataSet.AddEntry(Event(eventObs));
    }
    
    // name the observables for the save
    std::vector<std::string> names;
    names.push_back("var0");
    names.push_back("var1");
    names.push_back("var2");    
    names.push_back("var3");
    origDataSet.SetObservableNames(names);
    
    IO::SaveDataSet(origDataSet, "data_set_io_root_test.h5");
    OXSXDataSet* loadedSet = IO::LoadDataSet("data_set_io_root_test.h5");
    size_t nEntries = origDataSet.GetNEntries();

    REQUIRE(origDataSet.GetObservableNames() == loadedSet->GetObservableNames());


    
    REQUIRE(origDataSet.GetNEntries() == loadedSet->GetNEntries());
    REQUIRE(origDataSet.GetNObservables() == loadedSet->GetNObservables());

	  REQUIRE(loadedSet->GetNEntries() == nEntries);
    

    REQUIRE(origDataSet.GetEntry(0).GetData() == loadedSet->GetEntry(0).GetData());
    REQUIRE(origDataSet.GetEntry(nEntries -1).GetData() == loadedSet->GetEntry(nEntries -1).GetData());
    
    
    remove("data_set_io_root_test.h5");
}
