/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.entity;

import altr.strategies.Strategy;
import java.util.concurrent.atomic.AtomicInteger;

/**
 *
 * @author Vitaly
 */
public class Person implements Cloneable {
    
    private static AtomicInteger uniqueId = new AtomicInteger();
    
    private long id;
    private double money;
    private double initialMoney;
    private boolean active;
    private double ego;
    private double group;
    private double altr;
    private String name;
    private Strategy strategy;

    public long getId() {
        return id;
    }    

    public void setId(long id) {
        this.id = id;
    }
    
    public double getMoney() {
        return money;
    }

    public void setMoney(double money) {
        this.money = money;
    }

    public boolean isActive() {
        return active;
    }

    public void setActive(boolean active) {
        this.active = active;
    }

    public double getEgo() {
        return ego;
    }

    public void setEgo(double ego) {
        this.ego = ego;
    }

    public double getGroup() {
        return group;
    }

    public void setGroup(double group) {
        this.group = group;
    }

    public double getAltr() {
        return altr;
    }

    public void setAltr(double altr) {
        this.altr = altr;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Strategy getStrategy() {
        return strategy;
    }

    public void setStrategy(Strategy strategy) {
        this.strategy = strategy;
    }

    public double getInitialMoney() {
        return initialMoney;
    }
    
    @Override
    public Object clone() throws CloneNotSupportedException {
        Person p = (Person) super.clone(); 
        p.active = active;
        p.altr = altr;
        p.ego = ego;
        p.group = group;
        p.id = uniqueId.getAndIncrement();;
        p.money = money;
        p.name = name;
        return p;
    }

    @Override
    public String toString() {
        //return "Person{" + "id=" + id + ", money=" + money + ", active=" + active + ", ego=" + ego + ", group=" + group + ", altr=" + altr + ", name=" + name + '}';
        return "id=" + id + " money=" + money;
    }

    public Person(double money, boolean active, double ego, double group, double altr, String name, Strategy strategy) {
        this.id = uniqueId.getAndIncrement();
        this.initialMoney = money;
        this.money = money;
        this.active = active;
        this.ego = ego;
        this.group = group;
        this.altr = altr;
        this.name = name;
        this.strategy = strategy;
    }  
    
}
