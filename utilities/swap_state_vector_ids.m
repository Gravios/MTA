function stateVector = swap_state_vector_ids(stateVector,ind1,ind2);
tempa = stateVector==ind1;
tempb = stateVector==ind2;
stateVector(tempa) = ind2;
stateVector(tempb) = ind1;
