"""
for debugging the DRT_DP fitting method
"""


def local_standard_input(all_test_data):
    #    N_freqs, Z_exp = add
    N_freqs, Z_exp = choose_test(all_test_data, reduce=False)
    Z_exact, gamma_exact = standard_DRT(N_freqs)


#    print('distance_opt= ', distance_vec[index_opt])


#    analyze_results()


def train_model():
    model = vanilla_model()

    # initialize following variables
    zeta = torch.randn(N, N_zeta)
    loss_vec = np.array([])
    distance_vec = np.array([])
    lambda_vec = np.array([])

    # optimize the neural network
    learning_rate = 1e-5
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    _lambda = 100
    # max iterations
    max_iters = int(100001)
    gamma_NN_store = torch.zeros((max_iters, N_freqs))
    R_inf_NN_store = torch.zeros((max_iters, 1))
    C_const_NN_store = torch.zeros((max_iters, 1))

    for t in range(max_iters):
        # Forward pass: compute predicted y by passing x to the model.
        gamma = model(zeta)

        # Compute the loss
        loss = loss_fn(
            gamma, Z_exp_re_torch, Z_exp_im_torch, A_re_torch, A_im_torch, _lambda
        )
        # save it
        loss_vec = np.append(loss_vec, loss.item())

        # store gamma
        gamma_NN = gamma[:, 0:-1].detach().reshape(-1)
        gamma_NN_store[t, :] = gamma_NN

        # store R_inf
        R_inf_NN_store[t, :] = gamma[:, -1].detach().reshape(-1)

        # Compute the distance
        distance = math.sqrt(torch.sum((gamma_NN - gamma_exact_torch) ** 2).item())
        # save it
        distance_vec = np.append(distance_vec, distance)

        # and print it
        if not t % 10000:
            print("iter=", t, "; loss=", loss.item(), "; distance_NOT USED=", distance)

        # zero all gradients (purge any cache)
        optimizer.zero_grad()

        # compute the gradient of the loss with respect to model parameters
        loss.backward()

        # Update the optimizer
        optimizer.step()

    return distance_vec, loss_vec, gamma_NN_store, R_inf_NN_store


def analyze_results():
    index_opt = np.argmin(distance_vec)
    index_early_stop = np.flatnonzero(np.abs(np.diff(loss_vec)) < 5e-2)

    gamma_DIP_torch_opt = gamma_NN_store[index_opt, :]
    R_inf_DIP_torch_opt = R_inf_NN_store[index_opt, :]

    gamma_DIP_opt = gamma_DIP_torch_opt.detach().numpy()
    R_DIP_opt = R_inf_DIP_torch_opt.detach().numpy()

    if len(index_early_stop):
        gamma_DIP_torch_early_stop = gamma_NN_store[index_early_stop[0], :]
        gamma_DIP = gamma_DIP_torch_early_stop.detach().numpy()
        R_DIP = R_inf_NN_store[index_early_stop[0], :]
        R_DIP = R_DIP.detach().numpy()
        print("distance_early_stop = ", distance_vec[index_early_stop[0]])
    else:
        gamma_DIP = gamma_DIP_opt
        R_DIP = R_DIP_opt

    Z_DIP = R_DIP + np.matmul(A_re, gamma_DIP) + 1j * np.matmul(A_im, gamma_DIP)

    print("total number parameters = ", compute_DRT.count_parameters(model))

    print("distance_opt= ", distance_vec[index_opt])

    plot_loss()
    plot_err_iter()
    plot_imp()
    plot_DRT()
