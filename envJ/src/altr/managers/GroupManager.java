/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.managers;

import altr.entity.Group;
import altr.entity.Membership;
import altr.entity.Person;
import java.util.ArrayList;
import java.util.Collection;
import altr.Environment;

/**
 *
 * @author Vitaly
 */
public class GroupManager {
    private final Environment env;

    public GroupManager(Environment env) {
        this.env = env;
    }

    public Collection<Group> getGroups() {
        return env.getGroups();
    }

    public void setGroups(Collection<Group> groups) {
        env.setGroups(groups);
    }

    public Collection<Membership> getMemberships() {
        return env.getMemberships();
    }

    public void setMemberships(Collection<Membership> memberships) {
        env.setMemberships(memberships);
    }   
    
    public void addPeopleToGroup(Collection<Person> people, long groupId) {
        for(Person p: people) {
            getMemberships().add(new Membership(p.getId(), groupId));
        }        
    }   
    
    public void addGroup(Group g) {
        getGroups().add(g);
    }
     
    public Collection<Long> getMembersByGroupId(long groupId){
        if (getGroups() == null || getMemberships() == null) return null;
        
        Collection<Long> members = new ArrayList();
        for(Membership m: getMemberships()) {
            if (m.getGroupId() == groupId) members.add(m.getPersonId());
        }
        
        return members;
    }
    
    public Collection<Long> getGroupsByPersonId(long personId){
        if (getGroups() == null || getMemberships() == null) return null;
        
        Collection<Long> personGroups = new ArrayList();
        for(Membership m: getMemberships()) {
            if (m.getPersonId() == personId) personGroups.add(m.getGroupId());
        }
        
        return personGroups;
    }   
    
    public void clear() {
        getGroups().clear();
        getMemberships().clear();
    }
    
}