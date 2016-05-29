/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.strategies;

import altr.entity.Offer;
import altr.entity.Person;
import altr.managers.Analyzer;
import altr.managers.GroupManager;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Vitaly
 */
public class SimpleGroupStrategy implements Strategy {
    
    protected final GroupManager gM;
    private final double threshold;
    private final Collection<Long> acceptedGroupIds = new ArrayList<>();
    private final Collection<Long> declinedGroupIds = new ArrayList<>();
    private int votes = 0;
    private int maxVotesCount;

    public SimpleGroupStrategy(GroupManager gM) {
        this.gM = gM;
        this.threshold = 0;        
    }
    
    public SimpleGroupStrategy(GroupManager gM, double threshold) {
        this.gM = gM;
        this.threshold = threshold;        
    }

    @Override
    public boolean vote(Collection<Offer> offers, Collection<Person> people, long personId) {
        Long groupId = gM.getGroupsByPersonId(personId).iterator().next();
        //System.out.println(votes);
        if (acceptedGroupIds.contains(groupId)) {
            votes++;
            return true;
        } else if (declinedGroupIds.contains(groupId)) {
            votes++;
            return false;
        }
        Double groupAvg;
        try {
            groupAvg = Analyzer.getAverageOffer(offers, gM.getMembersByGroupId(groupId));
            boolean isAccepted = groupAvg > threshold;
            if (isAccepted) {
                acceptedGroupIds.add(groupId);
            } else {
                declinedGroupIds.add(groupId);
            }
            votes++;
            return isAccepted;
        } catch (Exception ex) {
            Logger.getLogger(SimpleGroupStrategy.class.getName()).log(Level.SEVERE, null, ex);
            return false;
        }        
    }
    
    public void reset() {
        votes = 0;
        acceptedGroupIds.clear();
        declinedGroupIds.clear();
    }
    
    public void init(){
        //reset();
        maxVotesCount = gM.getMemberships().size();
        //System.out.println(maxVotesCount);
    }
}
