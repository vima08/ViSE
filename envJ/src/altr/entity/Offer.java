/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.entity;

/**
 *
 * @author Vitaly
 */
public class Offer {
    private double ammount; 
    private long personId;

    public double getAmmount() {
        return ammount;
    }

    public void setAmmount(double ammount) {
        this.ammount = ammount;
    }

    public long getPersonId() {
        return personId;
    }

    public void setPersonId(long personId) {
        this.personId = personId;
    }

    public Offer(double ammount, long personId) {
        this.ammount = ammount;
        this.personId = personId;
    }    

    @Override
    public String toString() {
        return "Offer{" + "ammount=" + ammount + ", personId=" + personId + '}';
    }    
    
}
