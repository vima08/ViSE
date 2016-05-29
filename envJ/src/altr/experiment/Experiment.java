/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.experiment;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;

/**
 *
 * @author Vitaly
 */
public class Experiment {
    private long id;
    private int iterationNumber;
    private String description;
    private Date date;
    private String name;
    private Collection<Stage> stages;        

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

    public Date getDate() {
        return date;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Collection<Stage> getStages() {
        return stages;
    }

    public void setStages(Collection<Stage> stages) {
        this.stages = stages;
    }

    public int getIterationNumber() {
        return iterationNumber;
    }

    public void setIterationNumber(int iterationCount) {
        this.iterationNumber = iterationCount;
    }    
    
    public int getStepNumber() {
        int stepNumber = 0;
        for(Stage s: getStages()) {
            stepNumber += s.getStepNum();
        }
        return stepNumber;
    }

    public Experiment(long id, int iterationCount, String description, String name) {
        this.id = id;
        this.iterationNumber = iterationCount; 
        this.description = description;
        this.date = new Date();
        this.name = name;
        this.stages = new ArrayList();
    }  
    
    public void addStage(Stage stage) {
        stages.add(stage);
    }

    @Override
    public String toString() {
        return "Experiment{" + "id=" + id + ", description=" + description + ", date=" + date + ", name=" + name + ", stages=" + stages + '}';
    }
    
}
