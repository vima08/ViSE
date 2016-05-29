/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.strategies;

import altr.entity.Offer;
import altr.entity.Person;
import java.util.Collection;

/**
 *
 * @author Vitaly
 */
public class SimpleEgoisticStrategy implements Strategy {
    
    private final double threshold;
    
    public SimpleEgoisticStrategy() {
        this.threshold = 0;
    }

    public SimpleEgoisticStrategy(double threshold) {
        this.threshold = threshold;
    }

    @Override
    public boolean vote(Collection<Offer> offers, Collection<Person> people, long personId) {
        for(Offer offer: offers) {
            if (offer.getPersonId() == personId) {
                return (offer.getAmmount() > threshold);
            }
        }
        return false;
    }    
}
