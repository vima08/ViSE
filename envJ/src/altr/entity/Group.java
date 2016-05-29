/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.entity;

import java.util.concurrent.atomic.AtomicInteger;

/**
 *
 * @author Vitaly
 */
public class Group {
    private static AtomicInteger uniqueId = new AtomicInteger();
    
    private long id;
    private String description;
    private String name;

    public long getId() {
        return id;
    }

    public void setId(long id) {
        this.id = id;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Group(String description, String name) {
        this.id = uniqueId.getAndIncrement();
        this.description = description;
        this.name = name;
    }    
    
}
