#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mt19937-64.h"
#include "xoshiro256plus.h"
#include "splitmix64.h"

typedef struct Logger { // holds samples and other data 
	double** samples;
	int num_samples;
	int number_of_agents;
	int logging;
	int backup;
	uint32_t seed;
	unsigned long long max_iterations;
} Logger;

typedef struct Node { // these are linked up into a circular ring, holding all fitness values
	struct Node* left;      // left and right are in the ring neighbor ordering
	struct Node* right;
    struct Node* next;      // next and prev are in the ranks ordering -- will be NULL if not in shortlist
    struct Node* prev;
	// int rank; // let this be -1 if it's not in the shortlist
	double value;
} Node;

typedef struct Master { // contains pointers to min node, the threshold node of the shortlist, and some other data
    struct Node* min_node;  // careful with this -- it will hold the bottom of the shortlist, but that may not be the correct minimum if you remove stuff without inserting properly
    struct Node* threshold_node; // sim.
	double threshold;
    int short_k;
    int short_length;
    unsigned long long relisted;
    int number_of_nodes;
	unsigned long long* rank_replacement_frequencies;
    char* rng;
} Master;

double next_real(void){
    return (next() >> 11) * (1.0/9007199254740991.0);
}

double* create_rand_array(char* rng, int number_of_agents){ // given a seed and a size, creates array of random fitnesses
	// initialize the array of agents, stop if the malloc went wrong
	double* array = malloc(number_of_agents * sizeof(double));
	if (array == NULL){
		free(array);
		return NULL;
	}

	// populate fitness array with random values in [0,1)

    if (strcmp(rng, "mt") == 0){
        for(int i = 0; i < number_of_agents; i++){
            array[i] = genrand64_real1();
        }
    } else if (strcmp(rng, "xoshiro") == 0){
        for(int i = 0; i < number_of_agents; i++){
            array[i] = next_real();
        }
    }

	// and return
	return array;
}

Node* create_node(Node* left, Node* right, Node* prev, Node* next, double value){
    Node* node = malloc(sizeof(Node));
    if (node == NULL){return NULL;}
    node->left = left;
    node->right = right;
    node->prev = prev;
    node->next = next;
    node->value = value;
	//node->rank = rank;
    return node;
}

Node* get_min_ptr(Node* head){ // returns pointer to node with minimum value in ring
    Node* min_ptr = head;
    Node* next_node = head->right;
	while(next_node != head){
		if (next_node->value < min_ptr->value){
			min_ptr = next_node;
		}
        next_node = next_node->right;
	}
	return min_ptr;
}

Node* create_ring(double* array, int number_of_agents) { // given an array of values creates a circular linked list of nodes and returns pointer to min. 
    // create head node
	Node* head = create_node(NULL, NULL, NULL, NULL, array[0]);
	if (head == NULL){return NULL;}
	head->value = array[0];

	// create all nodes in the masterlist starting from the head and linking appropriately
	Node* prev_node = head;
	for (int i = 1; i < number_of_agents - 1; i++){
		Node* next_node = create_node(prev_node, NULL, NULL, NULL, array[i]);
		if (next_node == NULL){return NULL;}
		prev_node->right = next_node;
		prev_node = next_node;
	}
	Node* next_node = create_node(prev_node, head, NULL, NULL, array[number_of_agents - 1]);
	if (next_node == NULL){return NULL;}
	prev_node->right = next_node;
	head->left = next_node;

	return get_min_ptr(head);
}

void delete_ring(Node* head){ // frees all the nodes in the circular linked master list

	// first decouple head from its left neighbor
	head->left->right = NULL;
	head->left = NULL;
	
	// now walk along and free as we go
	Node* temp;
	while (head != NULL){
		temp = head;
		head = head->right;
		free(temp);
	}
}

Master* create_master(int number_of_nodes, int short_k, char* rng){
	Master* master = malloc(sizeof(Master));
	double* init_vals = create_rand_array(rng, number_of_nodes);
	master->min_node = create_ring(init_vals, number_of_nodes);;
	master->threshold_node = NULL; 									// fix this later
	master->threshold = 0;
	master->short_k = short_k;
    master->short_length = 0;
    master->relisted = 0;
	master->number_of_nodes = number_of_nodes;
	master->rank_replacement_frequencies = calloc(number_of_nodes, sizeof(unsigned long long));
    master->rng = rng;
	//for (int i=0;i<5*short_k;i++){master->rank_replacement_frequencies[i] = 0;}
	return master;
}

void delete_master(Master* master){
	delete_ring(master->min_node);
	//free(master->rank_replacement_frequencies);
	free(master);
}

void* ring_to_array(Node* head, double* array){ // given pointers to circular linked list and an array, fills array with ring values
	int i = 0;
	Node* current = head;
	array[i] = current->value;
	current = current->right;
	i++;
	while(current != head){
		array[i] = current->value;
		current = current->right;
		i++;
	}
}

int compare_nodes(const void* p, const void* q ){ // given pointers to nodes p and q, return -1 if p < q, reuturn 1 if p > q (in value), 0 otherwise. used in qsort
	int out = 0;
	double pval = (*(Node**)p)->value;
	double qval = (*(Node**)q)->value;
	if (pval < qval){
		out = -1;
	}
	if (pval > qval){
		out = 1;
	}
	return out;
}

void linkup_shortlist(Master* master){ // given pointer to master, link nodes via ->next and ->prev into ranks

	// create array of pointers to nodes in the ring
	Node** node_pointers = malloc(master->number_of_nodes*sizeof(Node*));

	int i = 0;
	Node* current = master->min_node;
	node_pointers[i] = current;
    current->prev = NULL;  // wipe any shortlist connections
    current->next = NULL;
	current = current->right;
	i++;
	while(current != master->min_node){ // until we wrap around to the beginning, write down all node pointers in turn
		node_pointers[i] = current;
        current->prev = NULL;  // wipe any shortlist connections
        current->next = NULL;
		current = current->right;
		i++;
	}

	// sort it
	qsort((void*)node_pointers, master->number_of_nodes, sizeof(Node*), compare_nodes); 

    // take now only the smallest
    current = node_pointers[0]; // pointer to minimum node
    current->next = node_pointers[1];
    current = current->next;
    for(int j = 1; j < 5*master->short_k - 1; j++){
        current->prev = node_pointers[j - 1];
        current->next = node_pointers[j + 1];
        current = current->next;
    }
    current->prev = node_pointers[5*master->short_k - 2]; // current is threshold now

    master->min_node = node_pointers[0];
    master->threshold_node = current;
	master->threshold = master->threshold_node->value;

	// don't need this now
	free(node_pointers);
    master->short_length = 5*master->short_k;
	


}

void link(Node* p, Node* q){ // sets p and q to point at each other in order p -> q in shortlist
    p->next = q;
    q->prev = p;
}

void step_down(double pval, Node** qref){ // given a pointer to q=p, steps q down in ranks LL until q < pval < q->next
	while (pval < (*qref)->prev->value){
		*qref = (*qref)->prev;
	}
	*qref = (*qref)->prev;
	//printf("%d\n", *qref);
}

void step_up(double pval, Node** qref){ // given a pointer to q=p, steps q up in ranks LL until q < pval < q->next
	while (pval > (*qref)->next->value){
		*qref = (*qref)->next;
	}
	//printf("%d\n", *qref);
}

void remove_from_shortlist(Node* to_remove, Master* master){ // takes a pointer to a node in the shortlist we want to remove (relink on either side)
    if (to_remove->prev != NULL || to_remove->next != NULL){ // don't do anything if not actually in shortlist
        if (to_remove->prev == NULL){ // i.e. I am the minimum
            master->min_node = to_remove->next; // second rank is new minimum now
            to_remove->next->prev = NULL; // second rank now has no prev
            to_remove->next = NULL; // and I have no next
        } else if (to_remove->next == NULL) { // i.e. I am the threshold
            master->threshold_node = to_remove->prev; // second to last is now new threhsold 
            to_remove->prev->next = NULL; // second to last rank now has no next
            to_remove->prev = NULL; // and I have no prev
        } else { // I'm just in the middle of the ranks somewhere
            link(to_remove->prev, to_remove->next); // my neighbors should point to each other
            to_remove->next = NULL; // and I should point to no one
            to_remove->prev = NULL; 
        }
        master->short_length--; // size went down by one
		master->threshold = master->threshold_node->value;
    }
} 

void insert_bottom(Master* master, Node* p){ // p is out of order, but is less than threshold and is positioned:    p -> O -> O -> ...
	double pval = p->value;
	if (pval > p->next->value){ // otherwise do nothing, p is still the true minimum
		// in this case p->next is definitely the new min
		master->min_node = p->next;
		master->min_node->prev = NULL;

		// now we walk along until q < p < q->next. Note q->next is never NULL since p is assumed less than threshold
		Node* q = p->next;
		step_up(pval, &q); // q now points to node just below where p should be inserted (at worst q == threshold->prev)
		link(p, q->next);
		link(q, p);	
	}
}

void insert_middle(Master* master, Node* p){ // p is out of order, but is less than threshold and is positioned:   ... -> O -> p -> O -> ...
	double pval = p->value;
	if (pval < master->min_node->value){ // if p should be new min, just unlink it and attach to bottom
		link(p->prev, p->next); // my neighbors should point to each other
		link(p, master->min_node); // the old min and I should point to each other
		p->prev = NULL; // I am left end now
		master->min_node = p; // and I am new min now
	} else { 
		// here I can assume min < p < threshold
		if (p->prev->value > pval || pval > p->next->value){ // otherwise do nothing since p is in the right spot
			link(p->prev, p->next); // my neighbors should point to each other in any case

			Node* q = p; // set a temporary pointer q
			if (pval < p->prev->value){ // either p's rank should be decreased
				step_down(pval, &q); 
			} else if (pval > p->next->value){ // or p's rank should be increased
				step_up(pval, &q); 
			}
			//printf("%d\n", q);
			// q now points to node just below where p should be inserted (most extreme cases are q = min and q = threshold->prev)
			link(p, q->next);
			link(q, p);
		}
	}
}

void insert_top(Master* master, Node* p){ // p is out of order, but is less than original threshold and is positioned:   ... -> O -> p
	// careful here -- this is assuming that p was the threshold node, then its value was changed, but the new value is still
	// less than that original threshold value (otherwise we removed it already and didn't call this). So this function checks
	// if p is still largest (i.e. p->prev < p < old_threshold) in which case nothing happens (master still points to p which 
	// is the correct new smaller threshold node now). Otherwise, p->prev should be the new threshold in master, and we need to 
	// re-insert p in the right spot.
	double pval = p->value;
	if (pval < p->prev->value){
		// p->prev should be the new threshold
		master->threshold_node = p->prev;
		master->threshold_node->next = NULL;

		// now check if I'm the new min first
		if (pval < master->min_node->value){ 
			link(p, master->min_node); // the old min and I should point to each other
			p->prev = NULL; // I am left end now
			master->min_node = p; // and I am new min now
		} else { 
			// otherwise step down
			Node* q = p->prev; 
			step_down(pval, &q); 
			// q now points to node just below where p should be inserted (most extreme cases are q = min and q = threshold->prev)
			link(p, q->next);
			link(q, p);
		}
	}
}

void insert_free(Master* master, Node* p){ // p is not in the shortlist but should be
	double pval = p->value;
	if (pval < master->min_node->value){ 
		link(p, master->min_node); // the old min and I should point to each other
		master->min_node = p; // and I am new min now
	} else { 
		// otherwise step up from min
		Node* q = master->min_node; 
		step_up(pval, &q); 
		// q now points to node just below where p should be inserted (most extreme cases are q = min and q = threshold->prev)
		link(p, q->next);
		link(q, p);
	}
}

void shortlist_insert(Master* master, Node* p){ // given a node p that is less than threshold value, but unknown position in shortlist (or not at all), places p into shortlist correctly
	if (p->prev != NULL || p->next != NULL){ // i.e. p is in the shortlist already	
		if (p->prev == NULL){ // i.e. NULL -> p -> O -> ... so p is the bottom node
			insert_bottom(master, p);
		} else if (p->next == NULL) { // i.e. ... -> O -> p -> NULL so p is the top node
			insert_top(master, p);
		} else { // I'm just in the middle of the ranks somewhere ... -> O -> p -> O -> ...
			insert_middle(master, p);
		}
	} else { // in this case p is not connected to shortlist yet
		insert_free(master, p);
		master->short_length++; // shortlist got a new node
	}
	master->threshold = master->threshold_node->value;
}

void shortlist_update(Master* master){ // performs a single update using shortlist method
    // generate three new values for {left, middle, right}
    double newvals[3];

    if (strcmp(master->rng, "mt") == 0){
        newvals[0] = genrand64_real1();
        newvals[1] = genrand64_real1();
        newvals[2] = genrand64_real1();
    } else if (strcmp(master->rng, "xoshiro") == 0){
        newvals[0] = next_real();
        newvals[1] = next_real();
        newvals[2] = next_real();
    }

    // get pointers to the nodes that will change
    Node* middle = master->min_node;
    Node* left = middle->left;
    Node* right = middle->right;

	// FIRST RECORD RANKS TO BE REPLACED

	master->rank_replacement_frequencies[0]++; // middle is rank 0 always

	// determine the current ranks of left and right and record 

	// first check if left and right are even in the short list
	int needed = 0; // will hold wheter 0, 1, or 2 ranks will need to be found
	if (left->prev != NULL){ // left is not the min, so must have a prev if in shortlist
		needed++;
	}
	if (right->prev != NULL){ // sim for right
		needed++;
	}

	// now step through seeking needed ranks
	int i = 1;
	Node* current = middle->next; // start from the second rank
	int found = 0; // holds how many ranks have been found already
	while (found < needed){
		if (current == left){
			master->rank_replacement_frequencies[i]++;
			found++;
		}
		if (current == right){
			master->rank_replacement_frequencies[i]++;
			found++;
		}
		current = current->next;
		i++;
	}

	// NOW DO THE REPLACEMENTS


	middle->value = newvals[1]; 
	// for the middle node, determine if it should remain in the shortlist
	if (middle->value >= master->threshold ){ // in this case, just remove the node from the shortlist, set the new min node, and update value
		master->min_node = middle->next; // second rank is new minimum now
		master->min_node->prev = NULL; // new min now has no prev
		middle->next = NULL; // and I have no next
		master->short_length--; // shortlist got smaller
	} else {
		// in here I'm guaranteed that the new value is less than the threshold and also that 'middle' is the bottom of shortlist
		insert_bottom(master, middle);
		// shortlist length did not change
	}
	master->threshold = master->threshold_node->value;

	// for the left and right nodes, determine if they should be in shortlist and inssert/reposition if so

	left->value = newvals[0];

	//remove_from_shortlist(left, master);
	if (left->value >= master->threshold){
		remove_from_shortlist(left, master); // just get rid of it
	} else {
		shortlist_insert(master, left); // takes care of all cases assuming value is < threshold
	}

	right->value = newvals[2];

	//remove_from_shortlist(right, master);
	if (right->value >= master->threshold){
		remove_from_shortlist(right, master); // just get rid of it
	} else {
		shortlist_insert(master, right); // takes care of all cases assuming value is < threshold
	}

	while (master->short_length > master->short_k*5){ // if ever we surpass 5k, bump off top nodes 
		remove_from_shortlist(master->threshold_node, master);
	}
}

int shortlist_method(Master* master){
    if (master->short_length < 3*master->short_k){
        linkup_shortlist(master);
        master->relisted++;
    }

    shortlist_update(master);

    return 1;
}

Logger* create_logger(int number_of_agents, int number_of_samples, unsigned long long max_iterations, int logging){ // creates a logger which will hold the samples collected and some other data
	// create new empty logger, stop if the malloc went wrong
	Logger* logger = malloc(sizeof(Logger));
	if (logger == NULL){
		return NULL;
	}

	logger->num_samples = number_of_samples;
	logger->number_of_agents = number_of_agents;
	logger->max_iterations = max_iterations;

	// initialize array of arrays to hold samples
	logger->samples = calloc(logger->num_samples, sizeof(double*));
	if (logger->samples == NULL){
		free(logger);
		return NULL;
	}
	for(int i = 0; i < logger->num_samples; i++){
		logger->samples[i] = calloc(number_of_agents, sizeof(double));
	}
	logger->logging = logging;
	return logger;
}

void delete_logger(Logger* logger){ // frees everything
	if (logger != NULL){
		free(logger->samples);
		free(logger);
	}
}

int rings_equal(Node* p, Node* q){ // compares two rings to see if the fitness arrays match, assume same length
    if(p->value != q->value){return 0;}
	Node* next_p = p->right;
    Node* next_q = q->right;
	while(next_p != p){
		if(p->value != q->value){return 0;}
	    next_p = next_p->right;
        next_q = next_q->right;
	}
    return 1;
}

void pprint(double double_num){ // prints doubles with correct precision
	printf("%.*lf\n", 17, double_num);
}

void print_ring(Master* master){ // prints out all nodes in ring
    Node* head = master->min_node;
	pprint(head->value);
	Node* next_node = head->right;
	while(next_node != head){
		pprint(next_node->value);
		next_node = next_node->right;
	}
}

void print_shortlist(Master* master){ // prints out all shortlist nodes in rank order
    Node* head = master->min_node;
	pprint(head->value);
	Node* current = head->next;
	while(current != NULL){
		pprint(current->value);
		current = current->next;
	}
}

void print_array(double* array, int number_of_agents){
	for(int j = 0; j < number_of_agents; j++){
		pprint(array[j]);
	}
}

void print_rank_freq(Master* master){
	for(int j = 0; j < master->number_of_nodes; j++){
		printf("%llu\n",master->rank_replacement_frequencies[j]);
	}
}

int rank_lists_equal(Master* master1, Master* master2){ // given two masters, returns -n where n is the number of mismatches in the ranks list
	int out = 0;
	if(master1->number_of_nodes == master2->number_of_nodes){	
		for(int j = 0; j < master1->number_of_nodes; j++){
			if(master1->rank_replacement_frequencies[j] != master2->rank_replacement_frequencies[j]){
				out -= 1;
			}
		}
		return out;
	} else {
		return 1; // mismatched lengths 
	}
}

int run_shortlist(Master* master, Logger* logger){
	//init_genrand(run_seed);
	unsigned long long max_iterations = logger->max_iterations;
	int number_of_agents = logger->number_of_agents;
	int num_samples = logger->num_samples;

	clock_t start = clock(), diff;

	char path[300];
	double* snapshot_array = calloc(logger->number_of_agents, sizeof(double));
	uint32_t start_time = time(0);
	char* rng = master->rng ;
	uint32_t seed = logger->seed;
	
	unsigned long long i = 0; // current iteration
	unsigned long long k = 0; // sample counter
	while (i < max_iterations){
		i = i + shortlist_method(master);

		if (logger->backup && max_iterations > 100 && i % (max_iterations/20) == 0){
			// this is for taking a snapshot of the current information in case the simulation gets disrupted
			// currently does every 5% of the runtime. Want to fully reconstruct the sim from the snapshot if needed

			// probably should make this its own function?

			// make a filename for this iteration's data
			snprintf(path, 50, "/scratch/cafish/%ldbak_temp%d.txt", start_time,i*100/max_iterations);
			//snprintf(path, 300, "C:\\Users\\Cameron\\OneDrive\\Math\\Research\\Bak Sneppen\\Code\\C\\To Upload\\%ldbak_temp%d.txt", start_time,i*100/max_iterations);

			// save the current fitness values
			ring_to_array(master->min_node, snapshot_array);

			// open the file
			FILE *fp;
			fp = fopen(path, "a+");

			// current iteration
			fprintf(fp, "i = %llu\n", i);
			
			// write down seed and which rng we used
			if (strcmp(rng, "mt") == 0){
				fprintf(fp, "mt seed: %ld\n", seed);
			} else if (strcmp(rng, "xoshiro") == 0){
				fprintf(fp, "xoshiro seed: %ld\n", seed);
			}

			// record some info on the number of nodes, k, number of relists
			fprintf(fp, "n: %d\n", master->number_of_nodes);
			fprintf(fp, "k: %d\n", master->short_k);
			fprintf(fp, "times relisted: %llu\n", master->relisted);

			// record the current replacement frequencies list
			fprintf(fp, "first 3k replacement frequencies: ");
			fprintf(fp, "[");
			for(int j = 0; j < 3*master->short_k - 1; j++){
				fprintf(fp, "%llu,",master->rank_replacement_frequencies[j]);
			}
			fprintf(fp, "%llu",master->rank_replacement_frequencies[3*master->short_k - 1]);
			fprintf(fp, "]");
			fprintf(fp, "\n");


			// finally also write down the current fitnesses
			fprintf(fp, "snapshot: [");
			int j = 0;
			while(j < logger->number_of_agents - 1){
				fprintf(fp, "%.*lf", 17, snapshot_array[j]);
				fprintf(fp, "; ");
				j++;
			}
			fprintf(fp, "%.*lf", 17, snapshot_array[j]);
			fprintf(fp, "]");
			fclose(fp);
		}

		// after halfway, start recording samples (enough to have 500 total in the end)
		if (logger->logging == 1){
			if (i > max_iterations/2){
				if (i % (max_iterations/(2*num_samples)) == 0){
					//printf("sampled at i is %d, k is %d\n", i, k);
					ring_to_array(master->min_node, logger->samples[k]);
					k++;
				}									
			}
		}
	}

	diff = clock() - start;
	int msec = diff*1000 / CLOCKS_PER_SEC;
	return msec;
}

float estimate_threshold(Logger* logger, int bins, int more){ // bins should be a power of 10
    int nsamples = logger->num_samples;
    int pop = logger->number_of_agents;

    // first we want a histogram i.e. an array where the nth index holds a count of the values in the interval
    // of [n/b, (n+1)/b) for some number of bins b

    int* histogram = calloc(bins, sizeof(int)); // bins should be a power of 10
    double sample; // will hold current sample
    int placement;  // will hold the index sample [i,j] should go

    for (int i = 0; i < nsamples; i++){
        for (int j = 0; j < pop; j++){
            histogram[(int)(bins*logger->samples[i][j])]++;
            //printf("sample: %f\n", logger->samples[i][j]);
           // printf("placement: %d\n", (int)(bins*logger->samples[i][j]));
        }
    } 

    double* freqs = calloc(bins, sizeof(double));  // normalize histogram to have integral 1 (remember dx = 1/bins so we need to multiply by bins to get the correct height values)
    for(int j = 0; j < bins; j++){
		freqs[j] = (bins*(float)histogram[j])/(pop*nsamples);
	}


    // we estimate threshold1 by looking for the first x value with density above 1.5
    int threshold1 = 0;
    while(freqs[threshold1] < 1.5 && threshold1 <= bins){
        threshold1++;
    }
    printf("threshold1 (first value whose estimated density is above 1.5): %.*lf\n", (int)log10(bins), (float)threshold1/bins);

    // we make another estimate for the threshold by looking for the average height over x values in [.7,.9)
    float height = 0;
    for(int j = 7*bins/10; j < 9*bins/10; j++){
        height+= freqs[j];
    }
    height = height / (2*bins/10);
   // printf("avg height = %.*lf\n", (int)log10(bins), height);
    float threshold2 = (1 - (1/height));
    printf("threshold2 (via avg height in upper area): %.*lf\n", (int)log10(bins), threshold2);


    if (more){
        printf("sampled fitness distribution: [");
        for(int j = 0; j < bins-1; j++){
            printf("%f,", freqs[j]);
        }
        printf("%f]\n", freqs[bins]);
    }


    // printf("\n");
    // for(int j = 0; j < bins; j++){
	// 	printf("%d, ", histogram[j]);
	// }


    free(histogram);
    free(freqs);
}

void print_out(Master* master, Logger* logger, int runtime, int more){
    printf("n: %d\n", master->number_of_nodes);
    printf("k: %d\n", master->short_k);
    printf("seconds elapsed: %d\n", runtime/1000);
    printf("times relisted: %llu\n", master->relisted);
    printf("first 3k replacement frequencies: ");
    printf("[");
	for(int j = 0; j < 3*master->short_k - 1; j++){
		printf("%llu,",master->rank_replacement_frequencies[j]);
	}
	printf("%llu",master->rank_replacement_frequencies[3*master->short_k - 1]);
	printf("]");
	printf("\n");
    estimate_threshold(logger, 1000, more);
    if (more == 1){
        printf("samples = ");
        printf("[");
        int i = 0;
        while(i < logger->num_samples - 1){
            printf("\'[");
            int j = 0;
            while(j < logger->number_of_agents - 1){
                printf("%.*lf", 17, logger->samples[i][j]);
                printf("; ");
                j++;
            }
            printf("%.*lf", 17, logger->samples[i][j]);
            printf("]\',");
            i++;
        }
        printf("\'[");
        int j = 0;
        while(j < logger->number_of_agents - 1){
            printf("%.*lf", 17, logger->samples[i][j]);
            printf("; ");
            j++;
        }
        printf("%.*lf", 17, logger->samples[i][j]);
        printf("]\'");
        printf("]\n");
    }
}

int main(int argc, char** argv){ // args should be {population, k, num_samples, rng} where 'rng' is 'mt' or 'xoshiro'

    int number_of_agents;
    unsigned long long max_iterations;
    int short_k; 
    int number_of_samples;
    char* rng;

    int more = 1; // if 1 then full samples information is printed

    if (argc == 5){ // if we specified arguments then set them
        number_of_agents = atoi(argv[1]);
        max_iterations = 2*(pow(number_of_agents,3));
        short_k = atoi(argv[2]);
        number_of_samples = atoi(argv[3]);	
        if (strcmp(argv[4], "mt") == 0){
            rng = "mt";
        } else if (strcmp(argv[4], "xoshiro") == 0){
            rng = "xoshiro";
        }
    } else { // otherwise some default options
        number_of_agents = 100;
        max_iterations = 2*(pow(number_of_agents,3));
        short_k = 5;
        number_of_samples = 10;	
        rng = "mt";
    }

	Logger* logger = create_logger(number_of_agents, number_of_samples, max_iterations, 1);
	logger->backup = 1;

    // set up the rng that we chose (similar strcmps appear in the building and running as well)
    if (strcmp(rng, "mt") == 0){
        uint32_t seed = time(0);
        init_genrand64(seed);
        printf("mt seed: %ld\n",seed);
		logger->seed = seed;
        //pprint(genrand64_real1());
    } else if (strcmp(rng, "xoshiro") == 0){
        tosplit = time(0);
        printf("xoshiro seed: %ld\n", tosplit);
        s[0] = split();
        s[1] = split();
        s[2] = split();
        s[3] = split();
		logger->seed = tosplit;
        //pprint(next_real());
    }
     
    Master* master = create_master(number_of_agents, short_k, rng);
    linkup_shortlist(master);
    
	
    int time = run_shortlist(master, logger);
    print_out(master, logger, time, more);
    
    delete_master(master);
    delete_logger(logger);	
}

