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
public interface Strategy {
    
    public boolean vote(Collection<Offer> offers, Collection<Person> people, long personId);
}
