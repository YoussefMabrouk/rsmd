/************************************************
 *                                              *
 *                rs@md                         *
 *    (reactive steps @ molecular dynamics )    *
 *                                              *
 ************************************************/
/* 
 Copyright 2020 Myra Biedermann
 Licensed under the Apache License, Version 2.0 
*/

#include "simulatorRate.hpp"

//
// setup stuff specific to hybrid MC/MD simulation with rate 
// based acceptance criterion
//
void SimulatorRate::setup(const Parameters& parameters)
{
    rsmdLOG( "setting up the simulation world ..." );

    // setup general stuff
    SimulatorBase::setup(parameters);

    // setup specific stuff
    rsFrequency = parameters.getOption("reaction.frequency").as<REAL>();

    // check statistics file and write header
    STATISTICS_FILE << std::setw(10) << "# cycle"
                    << std::setw(15) << "# candidates"
                    << std::setw(15) << "# accepted"
                    << std::setw(15) << "# attempted" << '\n';

    rsmdLOG( "... setup done, time to start the simulation!" );
    rsmdLOG( std::flush << std::setprecision(3) );
}

// do reactive step
//
//
void SimulatorRate::reactiveStep()
{
    std::vector<int> nReactionsAttempted(8, 0);
    std::vector<int> nReactionsAccepted(8, 0);
    std::stringstream accepted_string;
    std::stringstream attempted_string;
    int ntotalaccepted = 0;
    int ntotalattempted = 0;
    int i;
    std::vector<ReactionCandidate> acceptedCandidates {};
    std::unordered_map<std::string, int> candidateTypes {};

    // search for candidates
    universe.update(lastReactiveCycle);
    auto candidates = universe.CellSearchReactionCandidates(); //searchReactionCandidates();  // returns shuffled vector of reaction candidates
    STATISTICS_FILE << std::setw(10) << currentCycle << std::setw(15) << candidates.size();
    if( candidates.size() > 0 )
    {
        rsmdLOG( "... found " << candidates.size() << " potential reaction candidates" );
        // go through candidates and react them if accepted
        for( auto& candidate: candidates )
        {
            if( universe.isAvailable(candidate) )
            {
                for (i=0; i<universe.getReactionTemplates().size(); i++)
                {
                    if (candidate.reaction_name()==universe.getReactionTemplates()[i].getName())
                    {
                        ++ nReactionsAttempted[i];
                    }
                }
                if( acceptance(candidate) )
                {
                    universe.react(candidate);
                    acceptedCandidates.push_back(candidate);
                    rsmdLOG( "... reacted candidate " << candidate.shortInfo() );
                    for (i=0; i<universe.getReactionTemplates().size(); i++)
                    {
                        if (candidate.reaction_name() == universe.getReactionTemplates()[i].getName())
                        {
                            ++ nReactionsAccepted[i];
                        }
                    }
                }
            }
            else
            {
                rsmdDEBUG( candidate.shortInfo() << " is no longer available for reaction" );
            }
            candidateTypes.try_emplace( candidate.getName(), 0 );
            candidateTypes[candidate.getName()] += 1;
        }     
        
        for(std::vector<int>::iterator it = nReactionsAccepted.begin(); it != nReactionsAccepted.end(); ++it)
        ntotalaccepted += *it;
        for(std::vector<int>::iterator it = nReactionsAttempted.begin(); it != nReactionsAttempted.end(); ++it)
        ntotalattempted += *it;
        
        std::copy(nReactionsAccepted.begin(), nReactionsAccepted.end(), std::ostream_iterator<int>(accepted_string, " "));  
        std::copy(nReactionsAttempted.begin(), nReactionsAttempted.end(), std::ostream_iterator<int>(attempted_string, " "));
        
        STATISTICS_FILE << std::setw(50) << accepted_string.str().c_str() << std::setw(50) << attempted_string.str().c_str();
        
        // relaxation
        if( ntotalaccepted > 0 )
        {
            universe.write(currentCycle);
            rsmdLOG( "... reacted " << ntotalaccepted << " out of " << ntotalattempted << " available candidates (out of " << candidates.size() << " candidates)" );
            rsmdLOG( "... candidates were: ")
            for( const auto& pair : candidateTypes) 
            {
                rsmdLOG( "... " << pair.second << " " << pair.first ); 
            }
            if( mdEngine->runRelaxation(currentCycle) )
            {
                rsmdLOG( "... relaxation succeeded!" );
                lastReactiveCycle = currentCycle;
                ++ nCyclesReaction;
                // read configuration after relaxation and check if sensible
                universe.readRelaxed(currentCycle);
                for(auto& accepted: acceptedCandidates)
                {
                    universe.checkMovement(accepted);
                }
            } 
            else
            {
                rsmdWARNING( "... relaxation failed, stepping out!" );
                raise(SIGABRT);
            }
        }
        else
        {
            rsmdLOG( "... no candidates were accepted" );
            ++ nCyclesNoReaction;
        }
        
    }
    else
    {
        rsmdLOG( "...found no candidates");
        ++ nCyclesNoReaction;
    }
    STATISTICS_FILE << '\n' << std::flush;
}


//
// check acceptance
//
bool SimulatorRate::acceptance(const ReactionCandidate& candidate)
{
    REAL random = enhance::random(0.0, 1.0);
    REAL condition = rsFrequency * candidate.getCurrentReactionRateValue(); 
    rsmdDEBUG( "checking acceptance for candidate " << candidate.shortInfo() );
    rsmdDEBUG( "condition = " << rsFrequency << "*" << candidate.getCurrentReactionRateValue() << "=" << condition);
    if( random < condition )
    {
        rsmdDEBUG( "candidate accepted: " << random << " < " << condition );
        return true;
    }
    else 
    {
        rsmdDEBUG( "candidate rejected: " << random << " !< " << condition );
        return false;
    }
}



//
// finish & clean up
//
void SimulatorRate::finish() 
{
    STATISTICS_FILE.close();

    rsmdLOG( "" );
    rsmdLOG( "finished rs@md simulation" );
    rsmdLOG( "total " << (nCyclesReaction + nCyclesNoReaction) << " cycles have been performed:" );
    rsmdLOG( "      " << nCyclesReaction << " with reactions" );
    rsmdLOG( "      " << nCyclesNoReaction << " without reaction" );
    rsmdLOG( "      " << nCyclesFailedFirstRelaxation << " failed during the first relaxation attempt" );
    rsmdLOG( "" << std::flush );
}

