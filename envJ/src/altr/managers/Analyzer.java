/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.managers;

import altr.entity.Offer;
import altr.entity.Person;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author Vitaly
 */
public class Analyzer {
    
    private static Comparator<Person> comparator = new Comparator<Person>() {
        @Override
        public int compare(Person p1, Person p2) {
            return (int) (p1.getMoney() - p2.getMoney()); 
        }
    };
    
    public static double getAverageMoney(Collection<Person> people) {
        double sum = 0;
        for(Person p: people) {
            sum += p.getMoney();
        }
        return sum / people.size();
    }   
    
    public static double getAverageMoney(Collection<Person> people, Collection<Long> personIds) throws Exception {
        int count = 0;
        double sum = 0;
        for(Person p: people) {
            if (personIds.contains(p.getId())) {
                count++;
                sum += p.getMoney();
            }
        }
        if (count == personIds.size()) {
            return sum / count;
        } else {
            throw new Exception("getAverageMoney: count != personIds.size()");
        }
    } 
    
    public static double getAverageOffer(Collection<Offer> offers, Collection<Long> personIds) throws Exception {
        int count = 0;
        double sum = 0;
        for(Offer o: offers) {
            if (personIds.contains(o.getPersonId())) {
                count++;
                sum += o.getAmmount();
            }
        }
        if (count == personIds.size()) {
            return sum / count;
        } else {
            throw new Exception("getAverageOffer: count != personIds.size()");
        }
    }

    public static double getAverageOffer(Collection<Offer> offers) {
        int count = 0;
        double sum = 0;
        for(Offer o: offers) {
            count++;
            sum += o.getAmmount();
        }
        return sum / count;
    }
    
    public static long getPositiveOffersCount(Collection<Offer> offers) {
        long count = 0;        
        for(Offer o: offers) {
            if (o.getAmmount() > 0) 
                count++;
        }
        return count;
    }
    
    public static double getPositiveOffersProportion(Collection<Offer> offers) {
        double count = 0;        
        for(Offer o: offers) {
            if (o.getAmmount() > 0) 
                count += 1;
        }
        return count/offers.size();
    }
    
    public static List<Person> sort(Collection<Person> people) {
        List<Person> rezult = new ArrayList(people);
        Collections.sort(rezult, comparator);
        return rezult;
    }
}
