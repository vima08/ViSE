/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.strategies;

import altr.entity.Offer;
import altr.entity.Person;
import altr.managers.Analyzer;
import java.util.Collection;

/**
 *
 * @author Vitaly
 */
public class SimpleAltruisticStrategy implements Strategy {
    
    //private final double leftBound;
    //private final double rightBound;
    
//    public SimpleAltruisticStrategy(double lb, double rb) {
//        this.leftBound = lb;
//        this.rightBound = rb;
//    }

    private final double threshold;
    private Boolean result = null;

    public SimpleAltruisticStrategy(double threshold) {
        this.threshold = threshold;
    }

    public SimpleAltruisticStrategy() {
         this.threshold = 0;
    }

    @Override
    public boolean vote(Collection<Offer> offers, Collection<Person> people, long personId) {
        if (result == null) {
            Double socialAvg = Analyzer.getAverageOffer(offers);
            result = (socialAvg > threshold);
        }
        return result;
    }    

    public void reset() {
        result = null;
    }
}
