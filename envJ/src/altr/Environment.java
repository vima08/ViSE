/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr;

import altr.entity.Group;
import altr.entity.Membership;
import altr.entity.Person;
import java.util.ArrayList;
import java.util.Collection;
import java.util.concurrent.atomic.AtomicInteger;

/**
 *
 * @author Vitaly
 */
public class Environment {

    private static AtomicInteger uniqueId = new AtomicInteger();

    private long id;
    private String name;
    private Collection<Group> groups;
    private Collection<Membership> memberships;
    private Collection<Person> people;

    public Environment(String name) {
        this.id = uniqueId.getAndIncrement();
        this.name = name;
        this.groups = new ArrayList<>();
        this.memberships = new ArrayList<>();
        this.people = new ArrayList<>();
    }

    public long getId() {
        return id;
    }

    public void setId(long id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Collection<Group> getGroups() {
        return groups;
    }

    public void setGroups(Collection<Group> groups) {
        this.groups = groups;
    }

    public Collection<Membership> getMemberships() {
        return memberships;
    }

    public void setMemberships(Collection<Membership> memberships) {
        this.memberships = memberships;
    }

    public Collection<Person> getPeople() {
        return people;
    }

    public void setPeople(Collection<Person> people) {
        this.people = people;
    }

}
