import risefl_interface
import os
import base64


# initialize parameters
num_clients = 3
max_malicious_clients = 1
dim = 2
num_blinds_per_group_element = 2
random_normal_bit_shifter = 24

# check criteria of a model update u:
#   ||u - center|| <= norm_bound
#
# Requirement: norm_bound > 0
#
norm_bound = 2
center = [10, 20]
check_param = risefl_interface.CheckParamFloat(risefl_interface.CHECK_TYPE_SPHERE)
check_param.sphere_param.bound = norm_bound
check_param.sphere_param.center = risefl_interface.VecFloat(center)

num_norm_bound_samples = 1000

# num_norm_bound_samples = 3000

# num_norm_bound_samples = 9000

weight_updates_collection = [[0, 0], [10, 21], [10.3, 20.4], [10.6, 19.9]]

# non-private protocols
protocol_type = risefl_interface.PROTOCOL_TYPE_NON_PRIV_INT         # weight updates rounded to integers
# protocol_type = risefl_interface.PROTOCOL_TYPE_NON_PRIV_FLOAT     # use float32 weight updates

# if use risefl_interface.PROTOCOL_TYPE_NON_PRIV_INT, rounded weight updates to integers of bit-length weight_bits
weight_bits = 24

inner_prod_bound_bits = weight_bits + random_normal_bit_shifter + 4
max_bound_sq_bits = 2 * inner_prod_bound_bits + 100

# initialize server
server = risefl_interface.ServerInterface(
    num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
    weight_bits, random_normal_bit_shifter,
    num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
    risefl_interface.CHECK_TYPE_SPHERE,
    False, protocol_type)

server.initialize_new_iteration(check_param)

print("server.dim = " + str(server.dim))
print("server.weight_bits = " + str(server.weight_bits))

# initialize clients
clients = []
for i in range(num_clients+1):
    client = risefl_interface.ClientInterface(
        num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
        weight_bits, random_normal_bit_shifter,
        num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
        risefl_interface.CHECK_TYPE_SPHERE,
        i, risefl_interface.VecSignPubKeys(), risefl_interface.SignPrvKey(),
        False, protocol_type)
    # print("client.client_id = " + str(client.client_id))
    clients.append(client)

print("clients[0].dim = " + str(clients[0].dim))
print("clients[0].client_id = " + str(clients[0].client_id))

client_ids = range(1, num_clients+1)
print(client_ids)

# step 1: clients send messages to server
for i in client_ids:
    # weights_i = weight_updates_collection[i]
    # weights = risefl_interface.VecFloat(len(weights_i))
    # for j in range(len(weights_i)):
    #     weights[j] = weights_i[j]
    weights = risefl_interface.VecFloat(weight_updates_collection[i])
    # print(weights)
    # print(type(weights))
    server.receive_1(clients[i].send_1(check_param,
                                       weights), i)

server.finish_iteration()
print("***** iteration finished *****")
# finish one iteration, aggregate sum is in server.final_update_float
#   aggregate average is in server.final_update_float_avg

col_sums = [sum(x) for x in zip(*weight_updates_collection)]

for j in range(dim):
    assert (abs(num_clients * server.final_update_float_avg[j] - col_sums[j]) < 1e-4)

print("test python interface success")