package uk.ac.liv.mzidlib.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 *
 * @author SPerkins
 */
public class Java7Stream<T> {
        private Collection<T> collection = null;
        private boolean used = false;
        public static <T> Java7Stream<T> stream(Collection<T> coll) {
            Java7Stream<T> stream = new Java7Stream<>(coll);
            return stream;
        }
        
        private Java7Stream(Collection<T> collection) {
            this.collection = collection;           
        }
        
        public Java7Stream<T> filter(Java7Predicate predicate) {
            if (used) {
                throw new RuntimeException("This stream has already been used.");
            }
            
            used = true;
            List<T> newList = new ArrayList<>();
            for (T object : collection) {
                if (predicate.test(object)) {
                    newList.add(object);
                }
            }
            
                       
            return new Java7Stream<>(newList);
        }       
        
        public T withSmallestMappedDouble(Java7DoubleMapper<T> mapper) {
            if (used) {
                throw new RuntimeException("This stream has already been used.");
            }
            
            T smallest = null;
            for (T object : collection) {
                if (smallest == null || mapper.map(object) < mapper.map(smallest)) {
                    smallest = object;
                }
            }
            
            return smallest;
        }
        
        public Java7DoubleStream mapToDouble(Java7DoubleMapper<T> mapper) {
            if (used) {
                throw new RuntimeException("This stream has already been used.");
            }
            
            List<Double> newDoubleList = new ArrayList<>();
            for (T object : collection) {
                newDoubleList.add(mapper.map(object));
            }
            
            return Java7DoubleStream.doubleStream(newDoubleList);
        }
        
        public <S> Java7Stream<S> map(Java7Mapper<T, S> mapper) {
            if (used) {
                throw new RuntimeException("This stream has already been used.");
            }
            
            List<S> newList = new ArrayList<>();
            for (T object : collection) {
                newList.add(mapper.map(object));
            }
            
            return new Java7Stream<>(newList);
        }
        
        public List<T> toList() {
            if (used) {
                throw new RuntimeException("This stream has already been used.");
            }
            
            return (List<T>)collection;
        }       
        
        public Java7Optional<T> any() {
            if (used) {
                throw new RuntimeException("This stream has already been used.");
            }
            
            T object = null;
            if (collection.size() > 0) {
                object = collection.iterator().next();
            }
            
            final T obj = object;
            
            Java7Optional<T> opt = new Java7Optional<T>() {                
                @Override
                public T get() {
                    return obj;
                }

                @Override
                public boolean isEmpty() {
                    return obj == null;
                }
            };
            
            return opt;            
        }
    }