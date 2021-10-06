#!/usr/bin/env python
#


################
# Help message #
################
def help_message():
    print('''
            Command line flags: 
                -coord 
                        -> dir of the coordinates files 
                -ref 
                        -> name of the reference list files 
                -key 
                        -> name of the AMOEBA key file 
                -outdir 
                        -> dir of the output and intermediate files 
                -outfile 
                        -> name of the output file 
                -tinker 
                        -> full path of the tinker analyze.x exe 
                -gtol 
                        -> gradient tolerance to converge the optimization 
                            default: 0.01 
                -delta 
                        -> step size used to evaluate the numerical gradient 
                            default: 0.00001 
                -energy 
                        -> energy type, its value can be: 
                                tot: total interaction energy 
                                vdw: vdw interaction energy 
                                els: electrostatic interaction energy 
                                pol: polarization interaction energy 
                            default: tot 
                -optimizer 
                        -> Integer: the optimizer will be used 
                                0 : fmin_bfgs   : Quasi-Newton method 
                                1 : fmin        : Nelder-Mead Simplex algorithm 
                                2 : fmin_cg     : Non-linear (Polak-Ribiere) conjugate gradient algorithm 
                                3 : fmin_powell : Powell's (modified) level set method 
                -bound 
                        -> String: defines the lower and upper bounds of each parameter in optimization 
                                e.g. if one has 3 parameters in optimization, the format is: 
                                -bound '0.0, 5.0, 0.0, 6.0, 1.0, 7.0'
                -np 
                        -> number of parallel processes 
                -debug 
                        -> debug mode                                                               ''')


########################
# Test if a str an int #
########################
def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


#################
# Get QM energy #
#################
def get_qm(my_list, local_d, local_m1, local_m2, local_e, local_w):
    # read a list of QM energy into an array from a file
    # Input: my_list
    # returns: input struture and eqm array

    tmp = [line.split() for line in open(my_list)]

    npair = len(tmp)

    for n in range(0, npair):
         pair = tmp[n][0]
 
         if len(tmp[n]) == 4 :
             local_d.append(tmp[n][0])
             local_m1.append(tmp[n][1])
             local_m2.append(tmp[n][2])
             local_e.append(tmp[n][3])
             local_w.append(str("1.0"))
         elif len(tmp[n]) >= 5 :
             local_d.append(tmp[n][0])
             local_m1.append(tmp[n][1])
             local_m2.append(tmp[n][2])
             local_e.append(tmp[n][3])
             local_w.append(tmp[n][4])
             if len(tmp[n]) > 5 :
                 print("WARNING: more than 5 columns found for ", pair, "!!!")
                 print("         reading the first 5 ONLY!!!")
                 f_out.write( "WARNING: more than 5 columns found for " + str(pair) + "!!!\n" )
                 f_out.write( "         reading the first 5 ONLY!!!\n" )
         else:
             print("ERROR: No QM energy and/or input structure info found for ", pair, "!!!")
             exit()

    #return dimer mono1 mono2 eqm
    return


#####################################
# Get MM energy                     #
# quicker but have less safty check #
#####################################
#def get_mm(in_xyz, my_key_file):
def get_mm(tmp_args):
    # Input: in_xyz, my_key_file
    # returns: my_emm array

    import subprocess
    #import subprocess
    in_xyz = tmp_args[0]
    my_key_file = tmp_args[1]

    #my_emm = []
    my_emm = 0.0

    my_op1 = ""
    my_op2 = ""
    my_op3 = ""

    my_call1 = ""
    my_call2 = ""
    my_call3 = ""

    #npair = len(in_xyz)
    npair = 1
    #print npair

    my_tmp_output_dimer = out_dir + "tmp_dimer.tinker.out"
    my_tmp_output_monoa = out_dir + "tmp_monoA.tinker.out"
    my_tmp_output_monob = out_dir + "tmp_monoB.tinker.out"

    # excute the MM calculation, and collect the energy into an array

    my_grep_ene = ""

    if ene_type.lower() == "tot" :
        my_grep_ene = " em | grep  'Total Potential Energy' | awk '{print $5}'"
    elif ene_type.lower() == "vdw" :
        my_grep_ene = " em | grep 'Van der Waals' | awk '{print $4}'"
    elif ene_type.lower() == "els" :
        my_grep_ene = " em | grep 'Atomic Multipoles' | awk '{print $3}'"
    elif ene_type.lower() == "pol" :
        my_grep_ene = " em | grep 'Polarization' | awk '{print $2}'"

    #print my_grep_ene

    for n in range(0, npair):
         # get the pair name
         #pair  = in_xyz[n][0]
         #mono1 = in_xyz[n][1]
         #mono2 = in_xyz[n][2]
         pair  = in_xyz[0]
         mono1 = in_xyz[1]
         mono2 = in_xyz[2]

         xyz_di = coord_dir + pair
         xyz_m1 = coord_dir + mono1
         xyz_m2 = coord_dir + mono2

         my_op1 = ""
         my_op2 = ""
         my_op3 = ""

         my_call1 = str(tinker_exe + " " + xyz_di + " -k " + my_key_file + my_grep_ene)
         my_call2 = str(tinker_exe + " " + xyz_m1 + " -k " + my_key_file + my_grep_ene)
         my_call3 = str(tinker_exe + " " + xyz_m2 + " -k " + my_key_file + my_grep_ene)

         #print " running pair: ", n, " dimer..."
         #print " command: ", my_call1
         my_status1, my_op1 = subprocess.getstatusoutput(my_call1)
#6         #print " running pair: ", n, " monoA..."
#6         #print " command: ", my_call2
#6         my_status2, my_op2 = commands.getstatusoutput(my_call2)
#6         #print " running pair: ", n, " monoB..."
#6         #print " command: ", my_call3
#6         my_status3, my_op3 = commands.getstatusoutput(my_call3)
#6         #print my_status1
#6         #print my_status2
#6         #print my_status3

         #print subprocess.call(str(tinker_exe) + str(xyz_di) + " -k " + str(my_key_file) + my_grep_ene, shell=True)
         #print subprocess.call(str(tinker_exe) + str(xyz_m1) + " -k " + str(my_key_file) + my_grep_ene, shell=True)
         #print subprocess.call(str(tinker_exe) + str(xyz_m2) + " -k " + str(my_key_file) + my_grep_ene, shell=True)

         if my_op1 != '' :
#6             # if the monomer is too small to have 1-4 or non-bonded terms, their monomer vdw, els, and pol energy will not be printed in tinker
#6             if my_op2 == '' :
#6                 my_op2 = "0.0"
#6             if my_op3 == '' :
#6                 my_op3 = "0.0"
#6             etmp = float(my_op1) - float(my_op2) - float(my_op3)
             etmp = float(my_op1)
             #my_emm.append(str(etmp))
             my_emm = etmp
             if (idebug == 1) : 
                 print("   processing pair : %-4s  %-40s || E%-4s(AB, A, B): %12.3f %12.3f %12.3f" % (str(n), pair, str(ene_type), float(my_op1), float(my_op2), float(my_op3)))
                 f_out.write( "   processing pair : %-4s  %-40s || E%-4s(AB, A, B): %12.3f %12.3f %12.3f" % (str(n), pair, str(ene_type), float(my_op1), float(my_op2), float(my_op3)) + "\n")
         else:
             print("ERROR: something's wrong with the MM energy calculation for ", pair, "!!!")
             print("NO valide output!!!")
             exit()

    return my_emm


#####################################
## Get MM energy 2                  #
## Slower but have more safty check #
#####################################
#def get_mm_2(in_xyz, my_key_file):
#    # Input: in_xyz, my_key_file
#    # returns: my_emm array
#
#    import commands
#    #import subprocess
#
#    my_emm = []
#
#    npair = len(in_xyz)
#    #print npair
#
#    my_tmp_output_dimer = out_dir + "tmp_dimer.tinker.out"
#    my_tmp_output_monoa = out_dir + "tmp_monoA.tinker.out"
#    my_tmp_output_monob = out_dir + "tmp_monoB.tinker.out"
#
#    # excute the MM calculation, and collect the energy into an array
#
#    my_grep_ene = ""
#    my_awk_ene  = ""
#
#    if ene_type.lower() == "tot" :
#        my_grep_ene = " grep 'Total Potential Energy' "
#        my_awk_ene  = " | awk '{print $5}'"
#    elif ene_type.lower() == "vdw" :
#        my_grep_ene = " grep 'Van der Waals' "
#        my_awk_ene  = " | awk '{print $4}'"
#    elif ene_type.lower() == "els" :
#        my_grep_ene = " grep 'Atomic Multipoles' "
#        my_awk_ene  = " | awk '{print $3}'"
#    elif ene_type.lower() == "pol" :
#        my_grep_ene = " grep 'Polarization' "
#        my_awk_ene  = " | awk '{print $2}'"
#
#    #print my_grep_ene + my_awk_ene
#
#    my_tinker_error = " grep -i 'Total potential Energy' "
#
#    for n in xrange(0, npair):
#         # get the pair name
#         pair  = in_xyz[n][0]
#         mono1 = in_xyz[n][1]
#         mono2 = in_xyz[n][2]
#
#         xyz_di = coord_dir + pair
#         xyz_m1 = coord_dir + mono1
#         xyz_m2 = coord_dir + mono2
#
#         # Call tinker
#         #print xyz_di
#
#         my_op1 = ""
#         my_op2 = ""
#         my_op3 = ""
#
#         #my_status1, my_op1 = commands.getstatusoutput(tinker_exe + xyz_di + " -k " + my_key_file + my_grep_ene)
#         #my_status2, my_op2 = commands.getstatusoutput(tinker_exe + xyz_m1 + " -k " + my_key_file + my_grep_ene)
#         #my_status3, my_op3 = commands.getstatusoutput(tinker_exe + xyz_m2 + " -k " + my_key_file + my_grep_ene)
#         my_status1, my_op1 = commands.getstatusoutput(tinker_exe + xyz_di + " -k " + my_key_file + " em > " + my_tmp_output_dimer)
#         my_status2, my_op2 = commands.getstatusoutput(tinker_exe + xyz_m1 + " -k " + my_key_file + " em > " + my_tmp_output_monoa)
#         my_status3, my_op3 = commands.getstatusoutput(tinker_exe + xyz_m2 + " -k " + my_key_file + " em > " + my_tmp_output_monob)
#         #print my_status1
#         #print my_status2
#         #print my_status3
#
#         #print subprocess.call(str(tinker_exe) + str(xyz_di) + " -k " + str(my_key_file) + my_grep_ene, shell=True)
#         #print subprocess.call(str(tinker_exe) + str(xyz_m1) + " -k " + str(my_key_file) + my_grep_ene, shell=True)
#         #print subprocess.call(str(tinker_exe) + str(xyz_m2) + " -k " + str(my_key_file) + my_grep_ene, shell=True)
#
#         my_op1 = ""
#         my_op2 = ""
#         my_op3 = ""
#
#         my_op1 = commands.getoutput(my_tinker_error + my_tmp_output_dimer)
#         my_op2 = commands.getoutput(my_tinker_error + my_tmp_output_monoa)
#         my_op3 = commands.getoutput(my_tinker_error + my_tmp_output_monob)
#
#         if my_op1 != '' and my_op2 != '' and my_op3 != '' :
#             my_ene1 = ""
#             my_ene2 = ""
#             my_ene3 = ""
#
#             my_ene1 = commands.getoutput(my_grep_ene + my_tmp_output_dimer + my_awk_ene)
#             my_ene2 = commands.getoutput(my_grep_ene + my_tmp_output_monoa + my_awk_ene)
#             my_ene3 = commands.getoutput(my_grep_ene + my_tmp_output_monob + my_awk_ene)
#
#             if my_ene1 == '' :
#                 my_ene1 = "0"
#             if my_ene2 == '' :
#                 my_ene2 = "0"
#             if my_ene3 == '' :
#                 my_ene3 = "0"
#
#             etmp = float(my_ene1) - float(my_ene2) - float(my_ene3)
#             my_emm.append(str(etmp))
#         else:
#             print "ERROR: something's wrong with the MM energy calculation for ", pair, "!!!"
#             print "NO valide output!!!"
#             exit()
#
#    return my_emm


#################
# Cost function #
#################
def cost_func(my_parm_opt_tmp, input_xyz_tmp, *args):

    from multiprocessing import Pool
    #import sys
    #import pp

    # call get_qm and get_mm to get the energy array

    # Globle: line_no_opt[]
    # Globle: parm_name[]
    # Globle: parm_type[]
    # Globle: line_opt[]

    global n_cost_call
    global nproc
    # Input: my_parm_opt_tmp: array
    # Input: input_xyz_tmp: 2D list
    # Input: my_key_file (*args)
    # returns: cost

    my_ene_call = False

    my_parm_bound = 0.0
    my_upper = 0.0
    my_lower = 0.0

    if str(args[0]) == "" :
        # if no input key file, we make one
        # Call to generate the new key file
        #print "empty filename of key file"
        #my_key_file = "/home/qwang/py_opt_test/tmp.key"
        my_key_file = out_dir + str("tmp.key")
        key_gen(my_key_file, line_opt, my_parm_opt_tmp, line_no_opt, True)
    else :
        my_key_file = str(args[0])

    if str(args[1]).lower() == "ene_call":
        my_ene_call = True
    else :
        my_ene_call = False

    #eqm = get_qm(ref_list)
    #print input_xyz_tmp

    ## set up ppservers
    ## we don't have remote servers
    #ppservers = ()

    #if len(sys.argv) > 1:
    #    ncpus = int(sys.argv[1])
    #    # Creates jobserver with ncpus workers
    #    job_server = pp.Server(ncpus, ppservers=ppservers)
    #else:
    #    # Creates jobserver with automatically detected number of workers
    #    job_server = pp.Server(ppservers=ppservers)
    #
    #print "Starting pp with", job_server.get_ncpus(), "workers"

    # multiprocessing
    pool = Pool(processes=nproc)
    tmp_args = [(input_xyz_tmp[x], my_key_file) for x in range(0, len(input_xyz_tmp))]
    emm = pool.map(get_mm, tmp_args)

    # serial
    #emm = get_mm(input_xyz_tmp, my_key_file)

    if len(eqm) != len(emm) :
        print("ERROR: length of EQM and EMM are not equal!!!")
        exit()

    cost = 0.0
    
    # see if boundary is applied to the parameters
    if if_bound == True :
        if len(bound_list) != len(my_parm_opt_tmp) :
            print("ERROR: length of parameter list and the boundary list are not equal!!!")
            print(my_parm_opt_tmp)
            #print len(my_parm_opt_tmp)
            exit()

        my_parm_bound = 0.0
        for n in range(0, len(my_parm_opt_tmp)) :
            my_upper = float(bound_list[n][1])
            my_lower = float(bound_list[n][0])
            if float(my_parm_opt_tmp[n]) > my_upper :
                my_parm_bound = my_parm_bound + 10.0*(abs(float(my_parm_opt_tmp[n]) - my_upper)**2.0)
            elif float(my_parm_opt_tmp[n]) < my_upper :
                my_parm_bound = my_parm_bound + 10.0*(abs(float(my_parm_opt_tmp[n]) - my_lower)**2.0)

        # Calculate the cost
        # This doesn't seem work good, need more work here...
        for n in range(0, len(eqm)):
            cost = cost + abs((float(eqm[n]) - float(emm[n])) * float(weight[n]))*(1.0 + my_parm_bound)

    else :

        # Calculate the cost without boundary information
        for n in range(0, len(eqm)):
            cost = cost + abs((float(eqm[n]) - float(emm[n])) * float(weight[n]))
            #print weight[n]

    if my_ene_call : 
       
        if idebug == 1:
            print("")
            print("OUT: Steps of cost call: ", n_cost_call)
            print("OUT: Current parm set:", my_parm_opt_tmp)
            print("OUT: The cost of current parm set: ", cost)
            print("")
            print("    The individual pair energy ")
            print("    -------------------------------------------------------------------------")
            print("    Name %32s Weight %4s QM %12s MM (%4s ) " % ("", "", "", ene_type))
            print("")
            for n in range(0, len(eqm)):
                print("    %-35s %8.3f %13.6f %18.6f" % (str(dimer[n]), float(weight[n]), float(eqm[n]), float(emm[n])))
            print("    -------------------------------------------------------------------------")
            print("OUT: End of %10s cost call: " % str(n_cost_call))
            print("")
            print("")
            print("")

        f_out.write( "\n" )
        f_out.write( "OUT: Steps of cost call: " + str(n_cost_call) + "\n" )
        f_out.write( "OUT: Current parm set: " + str(my_parm_opt_tmp) + "\n" )
        f_out.write( "OUT: The cost of current parm set: " + str(cost) + "\n" )
        f_out.write( "\n" )
        f_out.write( "    The individual pair energy " + "\n" )
        f_out.write( "    -------------------------------------------------------------------------" + "\n" )
        f_out.write( "    Name %32s Weight %4s QM %12s MM (%4s ) " % ("", "", "", ene_type) + "\n" )
        f_out.write( "" + "\n" )
        for n in range(0, len(eqm)):
            f_out.write( "    %-35s %8.3f %13.6f %18.6f" % (str(dimer[n]), float(weight[n]), float(eqm[n]), float(emm[n])) + "\n" )
        f_out.write( "    -------------------------------------------------------------------------" + "\n" )
        f_out.write( "OUT: End of %10s cost call: " % str(n_cost_call) + "\n" )
        f_out.write( "" + "\n" )
        f_out.write( "" + "\n" )
        f_out.write( "" + "\n" )

        n_cost_call = n_cost_call + 1

    return cost


######################
# Numerical Gradient #
######################
def num_grad(my_parm_opt, input_xyz, *args):
    # Globle: line_no_opt = []
    # Globle: parm_name = []
    # Globle: parm_type = []
    # Globle: line_opt[]

    # Input: my_parm_opt: array
    # Input: input_xyz: 2D list

    import numpy

    # step size in calculate numerical gradient
    #delta_parm = 0.00001
    #delta_parm = 0.0000001

    # Gradient is an array with my_parm_opt long
    grad = []
    grad_array = []

    # local containers
    tmp_parm = []

    # copy the original parms to a new parm_plus list
    #for i in xrange(0, len(my_parm_opt)):
    #    #parm_plus.append(my_parm_opt[i])
    #    #parm_minus.append(my_parm_opt[i])
    #    tmp_parm.append(my_parm_opt[i])

    #numpy.asarray(parm_plus)
    #numpy.asarray(parm_minus)
    #numpy.asarray(my_parm_opt, float)
    #numpy.asarray(tmp_parm, float)
    tmp_parm = numpy.zeros(len(my_parm_opt))

    #my_plus_key  = str("/home/qwang/py_opt_test/plus.key")
    #my_minus_key = str("/home/qwang/py_opt_test/minus.key")
    my_plus_key  = out_dir + str("plus.key")
    my_minus_key = out_dir + str("minus.key")

    my_tmp_max = 0.0
    my_tmp_min = 0.0
    my_tmp_scale = 0.0

    if idebug == 1 :
        print("Get numberical gradient ...")
    f_out.write( "Get numerical gradient ..." + "\n" )

    for n in range(0, len(my_parm_opt)):

        for i in range(len(tmp_parm)):
            tmp_parm[i] = float(my_parm_opt[i])

        #delta_parm = delta_parm * tmp_parm[n]

        #print "before plus", tmp_parm[n]
        tmp_parm[n]  = tmp_parm[n]  + delta_parm
        #print "after plus", tmp_parm[n]

        # Call to generate the plus key file
        key_gen(my_plus_key, line_opt, tmp_parm, line_no_opt, False)
        tmp_cost_plus = cost_func(tmp_parm, input_xyz, my_plus_key, "not_ene")
        
        for i in range(len(tmp_parm)):
            tmp_parm[i] = float(my_parm_opt[i])

        #print "before minus", tmp_parm[n]
        tmp_parm[n] = tmp_parm[n] - delta_parm
        #print "after minus", tmp_parm[n]

        # Call to generate the minus key file
        key_gen(my_minus_key, line_opt, tmp_parm, line_no_opt, False)
        tmp_cost_minus = cost_func(tmp_parm, input_xyz, my_minus_key, "not_ene")
        
        gtmp = str((float(tmp_cost_plus) - float(tmp_cost_minus))/(2.0*delta_parm))

        grad.append(str(gtmp))

        #print "grad", grad

    if idebug == 1 :
        print("Numerical grad: ", grad)
        print("")
    f_out.write( "Numerical grad: " + str(grad) + "\n" )
    f_out.write( "\n" )

    grad_array = numpy.zeros(len(my_parm_opt))
    for i in range(0, len(grad)):
        grad_array[i] = float(grad[i])

    return grad_array


###################
# Read input parm #
###################
def read_parm(my_file):

    import numpy

    all_lines = []
    new_lines = []

    need_opt = []
    need_opt_2 = []
    no_opt = []

    my_opt_list = []
    my_sp_position = []  # 2D list
    my_parm_position = []  # 2D list

    my_if_opt = []
    my_if_stype = []
    my_stype_index = []
    my_stype_position = []
    my_if_stype = []
    my_sp_index = []
    my_sp_position = []

    my_n_pound_opt = 0

    # read the parameter file
    #new_lines = [line.split() for line in open(my_file)]
    for line in open(my_file):
        #all_lines.append(line.split('\n')[0])
        all_lines.append(str(line))

    for n in range(0, len(all_lines)):
        new_lines.append(str(all_lines[n]).split())

    nlines = len(new_lines)

    for n in range(0, nlines):
        tmp_iopt = False
        for nn in range(0, len(new_lines[n])):
            tmp_string = new_lines[n][nn].lower()
            if tmp_string[:4] == "#opt" :
                tmp_iopt = True
                break

        if tmp_iopt == True:
            need_opt.append(new_lines[n])
        else:
            no_opt.append(all_lines[n])

    my_sp_position = [[0] for i in range(len(need_opt))]  # 2D list
    my_parm_position = [[0] for i in range(len(need_opt))]  # 2D list

    my_if_opt = [0 for i in range(len(need_opt))]
    my_if_stype = [0 for i in range(len(need_opt))]
    my_stype_index = [0 for i in range(len(need_opt))]
    my_stype_position = [0 for i in range(len(need_opt))]
    my_sp_index = [0 for i in range(len(need_opt))]
    my_sp_position = [0 for i in range(len(need_opt))]

    # second loop
    for n in range(0, len(need_opt)):
        for nn in range(0, len(need_opt[n])):
            tmp_string = need_opt[n][nn].lower()
            if tmp_string[:4] == "#opt" :

                my_if_opt[n] = 1
                tmp_iread = 0
                tmp_line_np = 0

                tmp_string = need_opt[n][nn].lower()
                # one parameter
                if tmp_string == "#opt1":
                    tmp_line_np = 1
                    tmp_iread = nn
                # two parameters
                elif tmp_string == "#opt2":
                    tmp_line_np = 2
                    tmp_iread = nn
                # three parameters
                elif tmp_string == "#opt3":
                    tmp_line_np = 3
                    tmp_iread = nn
                # no parameters
                elif tmp_string == "#opt0":
                    tmp_line_np = 0
                    tmp_iread = 0
                    break
                else:
                    print("unknown #opt flag!")
                    print("need to be #opt1, #opt2, #opt3 or #opt0")
                    print(need_opt[n])
                    exit()

                # read in the column number of the params
                if tmp_line_np == 1:
                    if is_int(need_opt[n][tmp_iread + 1]):
                        tmp_id_1 = int(need_opt[n][tmp_iread + 1])
                        if tmp_id_1 < tmp_iread + 1 :
                            my_opt_list.append(need_opt[n][tmp_id_1 - 1])
                            my_parm_position[n][0] = tmp_id_1 - 1
                            my_n_pound_opt = my_n_pound_opt + 1
                        else:
                            print("the integer is larger than the cloumn number of the parameters!")
                            print(need_opt[n])
                            exit()
                    else:
                        print("An integer is expected after #opt1!")
                        print(need_opt[n])
                        exit()
                elif tmp_line_np == 2:
                    if (is_int(need_opt[n][tmp_iread + 1]) and 
                        is_int(need_opt[n][tmp_iread + 2]) ):
                        tmp_id_1 = int(need_opt[n][tmp_iread + 1])
                        tmp_id_2 = int(need_opt[n][tmp_iread + 2])
                        if (tmp_id_1 < tmp_iread + 1 and
                            tmp_id_2 < tmp_iread + 1) :
                            my_opt_list.append(need_opt[n][tmp_id_1 - 1])
                            my_opt_list.append(need_opt[n][tmp_id_2 - 1])
                            my_parm_position[n][0] = tmp_id_1 - 1
                            my_parm_position[n].append(tmp_id_2 - 1)
                            my_n_pound_opt = my_n_pound_opt + 2
                        else:
                            print("the integer is larger than the cloumn number of the parameters!")
                            print(need_opt[n])
                            exit()
                    else:
                        print("Two integers are expected after #opt2!")
                        print(need_opt[n])
                        exit()
                elif tmp_line_np == 3:
                    if (is_int(need_opt[n][tmp_iread + 1]) and 
                        is_int(need_opt[n][tmp_iread + 2]) and
                        is_int(need_opt[n][tmp_iread + 3])) :
                        tmp_id_1 = int(need_opt[n][tmp_iread + 1])
                        tmp_id_2 = int(need_opt[n][tmp_iread + 2])
                        tmp_id_3 = int(need_opt[n][tmp_iread + 3])
                        if (tmp_id_1 < tmp_iread + 1 and
                            tmp_id_2 < tmp_iread + 1 and
                            tmp_id_3 < tmp_iread + 1) :
                            my_opt_list.append(need_opt[n][tmp_id_1 - 1])
                            my_opt_list.append(need_opt[n][tmp_id_2 - 1])
                            my_opt_list.append(need_opt[n][tmp_id_3 - 1])
                            my_parm_position[n][0] = tmp_id_1 - 1
                            my_parm_position[n].append(tmp_id_2 - 1)
                            my_parm_position[n].append(tmp_id_3 - 1)
                            my_n_pound_opt = my_n_pound_opt + 3
                        else:
                            print("the integer is larger than the cloumn number of the parameters!")
                            print(need_opt[n])
                            exit()
                    else:
                        print("Three integers are expected after #opt3!")
                        print(need_opt[n])
                        exit()
            # End of #opt

            elif tmp_string[:6] == "#stype" :
                my_if_stype[n] = 1
                if my_if_opt[n] == 1 :
                    tmp_iread = 0

                    tmp_string = need_opt[n][nn].lower()
                    if tmp_string == "#stype":
                        tmp_iread = nn
                        if len(need_opt[n]) >= tmp_iread + 2 :
                            if is_int(need_opt[n][tmp_iread + 1]) :
                                tmp_id_1 = int(need_opt[n][tmp_iread + 1])
                                my_stype_index[n] = tmp_id_1

                                if len(need_opt[n]) >= tmp_iread + 3 :
                                    # we have an extra index after "#stype X"
                                    if is_int(need_opt[n][tmp_iread + 2]) :
                                        tmp_id_2 = int(need_opt[n][tmp_iread + 2])
                                        tmp_id_max = my_parm_position[n][0]
                                        tmp_id_min = my_parm_position[n][0]
                                        for kk in range(0, len(my_parm_position[n])):
                                            if tmp_id_max < my_parm_position[n][kk] :
                                                tmp_id_max = my_parm_position[n][kk]
                                            if tmp_id_min > my_parm_position[n][kk] :
                                                tmp_id_min = my_parm_position[n][kk]

                                        if tmp_id_2 - 1 > tmp_id_max or tmp_id_2 - 1 < tmp_id_min :
                                            print("index out of range; should be <= ", tmp_id_max + 1, " and >= ", tmp_id_min + 1)
                                            print(need_opt[n])
                                            exit()
                                        else :
                                            for ii in range(0, len(my_stype_index)):
                                                if ii != n :
                                                    if my_stype_index[n] > 0 and my_stype_index[n] == my_stype_index[ii] :
                                                        # This have the same type already, we need to remove it from the opt list
                                                        tmp_num_count = 0
                                                        for jj in range(0, len(my_parm_position[n])):
                                                            tmp_num_count = tmp_num_count + 1
                                                            #print "my_parm_position[n] ", my_parm_position[n]
                                                            if tmp_id_2 - 1 == my_parm_position[n][jj] :
                                                                del my_opt_list[len(my_opt_list) - len(my_parm_position[n]) + tmp_num_count - 1]
                                                                my_stype_position[n] = tmp_id_2
                                                                break
                                    else: 
                                        print("An integer is expected after #stype x !")
                                        print(need_opt[n])
                                        exit()
                                else :
                                    for ii in range(0, len(my_stype_index)):
                                        if ii != n :
                                            if my_stype_index[n] > 0 and my_stype_index[n] == my_stype_index[ii] :
                                                # This have the same type already, we need to remove it from the opt list
                                                if tmp_line_np == 1 :
                                                    my_opt_list.pop()
                                                    break
                                                elif tmp_line_np == 2 :
                                                    my_opt_list.pop()
                                                    my_opt_list.pop()
                                                    break
                                                elif tmp_line_np == 3 :
                                                    my_opt_list.pop()
                                                    my_opt_list.pop()
                                                    my_opt_list.pop()
                                                    break
                            else:
                                print("An integer is expected after #stype!")
                                print(need_opt[n])
                                exit()
                        else: 
                            print("An integer is expected after #stype !")
                            print(need_opt[n])
                            exit()
                else :
                    print("#stype must be accompanied by #opt !")
                    print(need_opt[n])
                    exit()
        # End of #stype

    return no_opt, need_opt, my_opt_list, my_parm_position, my_stype_position, my_stype_index, my_n_pound_opt
 


################################
# Generate additional key file #
################################
def key_gen(my_file, my_line_opt, p_opt, no_opt, my_ene_call):

    my_tmp_stype_holder = []

    my_tmp_stype_holder = [[0] for i in range(len(my_line_opt))]  # make it 2D list

    f = open(my_file, 'w')
    for n in range(0, len(no_opt)):
        #print no_opt[n]
        f.write(no_opt[n])

    np = 0
    for i in range(0, len(my_line_opt)):
        tmp_line = my_line_opt[i]
        #print tmp_line
        #print "parm position", parm_position[i]
        #print "sp position", sp_position[i]

        if stype_index[i] > 0 :
            tmp_same_type = False

            for jj in range(0, len(my_line_opt)) :
                if stype_index[i] == my_tmp_stype_holder[jj][0] :
                    tmp_same_type = True
                    break

            if tmp_same_type : 
                if stype_position[i] > 0 :
                    for ii in range(0, len(parm_position[jj])) :
                        if stype_position[i] - 1 == parm_position[jj][ii] :
                            tmp_index = parm_position[jj][ii]
                            tmp_line[tmp_index] = my_tmp_stype_holder[jj][ii + 1]
                        else:
                            tmp_index = parm_position[jj][ii]
                            #tmp_line[tmp_index] = my_tmp_stype_holder[jj][ii + 1]
                            tmp_line[tmp_index] = p_opt[np]
                            np = np + 1
                else:
                    for ii in range(0, len(parm_position[jj])) :
                        tmp_index = parm_position[jj][ii]
                        tmp_line[tmp_index] = my_tmp_stype_holder[jj][ii + 1]
            else :
                my_tmp_stype_holder[i][0] = stype_index[i]
                for ii in range(0, len(parm_position[i])) :
                    tmp_index = parm_position[i][ii]
                    tmp_line[tmp_index] = p_opt[np]
                    my_tmp_stype_holder[i].append(p_opt[np])
                    np = np + 1

        else:
            # no stype_index found
            for ii in range(0, len(parm_position[i])) :
                tmp_index = parm_position[i][ii]
                tmp_line[tmp_index] = p_opt[np]
                np = np + 1

        for ki in range(0, tmp_index+1):
            f.write('%-20s' % tmp_line[ki]),
            if my_ene_call :
                f_out.write('%-20s' % tmp_line[ki]),
            if idebug == 1 :
                print(('%-20s' % tmp_line[ki]), end=' ')

        for ki in range(tmp_index+1, len(tmp_line)-1):
            f.write('%-7s' % tmp_line[ki]),
            if my_ene_call :
                f_out.write('%-7s' % tmp_line[ki]),
            if idebug == 1 :
                print(('%-7s' % tmp_line[ki]), end=' ')

        f.write('%-7s\n' % tmp_line[len(tmp_line)-1])
        if my_ene_call :
            f_out.write('%-7s\n' % tmp_line[len(tmp_line)-1])
        if idebug == 1 :
            print(('%-7s' % tmp_line[len(tmp_line)-1]))

    f.close()
    return


################
# MAIN program #
################
#
#

#from scipy import optimize
import my_optimize
import numpy
import os
import sys
import time

# Globle
# Reference energy
dimer = []
mono1 = []
mono2 = []
eqm = []
weight = []
input_names = []

# lines that have "opt" sign
line_opt = []
line_no_opt = []

# parameters that need to be optimized
parm_name = []
parm_type = []
parm_opt_list = []
parm_opt = []
parm_position = [] # 2D list
stype_position = [] # 2D list
stype_index = [] # 1D list

# optimizer names
optimizer_name=['fmin_bfgs', 'fmin', 'fmin_cg', 'fmin_powell']

# Index for optimizer
# default : 0 : fmin_bfgs
optimizer = 0

# number of cost evaluation (except numerical gradient calls)
n_cost_call = 0

# number of processes
nproc = 1

# debug
idebug = 0

# number of parameters assigned the #opt
n_pound_opt = 0

# step size used to evaluate the numerical gradient
delta_parm = 0.00001

# if a boundary is applied to the parameters
bound_list = []  # 2D array

if_bound = False

tmp_bound_list = []



working_dir=os.getcwd()
working_dir=working_dir + str("/")

# Coordinate file path
# Default
coord_dir = "coord/"
parm_file = "force.key"
ref_file = "qm_list.dat"
out_dir = ""
out_file = "opt.log"

# TINKER exe
# Default
tinker_exe = "/home/pren/tinker6/source-2014-xz-node/analyze.x"

# Target energy type
#    tot - total interaction energy
#    vdw - vdW interaction energy
#    els - electrostatic interaction energy
#    pol - polarization interaction energy
# Default
ene_type = "tot"

# Gradient tolerance to converge the optimization
gtol = 0.01

xtol = 0.01
ftol = 0.01

# Timer
start_time = time.time()

# Read in command line arguments
if len(sys.argv) > 1 :
    for kk in range(0, len(sys.argv)) :
        if sys.argv[kk][:1] == "-" :
            if sys.argv[kk].lower() == "-coord" :
                coord_dir = str(sys.argv[kk + 1])
                if coord_dir[-1] != "/" :
                    # add the directory symbol to the variable
                    coord_dir = str(coord_dir) + "/"
            elif sys.argv[kk].lower() == "-ref" :
                ref_file = str(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-key" :
                parm_file = str(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-outdir" :
                if str(sys.argv[kk + 1]) == "./" :
                    out_dir = str(sys.argv[kk + 1])
                else :
                    out_dir = str(sys.argv[kk + 1]) + str("/")
            elif sys.argv[kk].lower() == "-outfile" :
                out_file = str(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-tinker" :
                tinker_exe = str(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-energy" :
                ene_type = str(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-debug" :
                idebug = 1
            elif sys.argv[kk].lower() == "-bound" :
                tmp_bound_list = str(sys.argv[kk + 1]).split(",")
            elif sys.argv[kk].lower() == "-optimizer" :
                optimizer = int(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-gtol" :
                gtol = float(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-delta" :
                delta_parm = float(sys.argv[kk + 1])
            elif sys.argv[kk].lower() == "-np" :
                nproc = int(sys.argv[kk + 1])
            #elif ( sys.argv[kk].lower() == "-h" or 
            #     sys.argv[kk].lower() == "--h" or
            #     sys.argv[kk].lower() == "-help" or
            #     sys.argv[kk].lower() == "--help" ):
            else:
                print("Unknown flags!")
                help_message()
                exit()

# Determine the energy type
if (ene_type.lower() != "tot" and ene_type.lower() != "vdw" and
    ene_type.lower() != "els" and ene_type.lower() != "pol") :
    print("Unknown -energy flag!", ene_type)
    help_message()
    exit()

if out_dir != "" :
    if out_dir[0] == "/" :
        # This is a full path
        out_dir = out_dir
    else :
        out_dir = working_dir + out_dir
else :
    out_dir = working_dir + out_dir
    
if coord_dir[0] != "/" :
    coord_dir = working_dir + coord_dir
parm_file = working_dir + parm_file

if ref_file[0] != "/" :
    ref_file = coord_dir + ref_file

out_file = out_dir + out_file

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

kk = 1
tmp_out_file = out_file
while os.path.isfile(out_file):
    kk = kk + 1
    out_file = tmp_out_file + "_" + str(kk)

# see if boundary is applied to the parameters
# transfer the 1D list to 2D
bound_list = numpy.array(tmp_bound_list).reshape(len(tmp_bound_list)/2,2)

if len(bound_list) > 0 :
    if_bound = True
    #print "bound applied"
    #print "bound: ", bound_list
    #print len(bound_list)
#else :
#    print "bound not applied"

# Open output file

f_out = open(out_file, 'w')

# Call to read the name of input strucutres for MM
# and the reference QM energy
get_qm(ref_file, dimer, mono1, mono2, eqm, weight)

#print dimer
#print mono1
#print mono2
#print eqm
#print weight

if len(dimer) == len(mono1) and len(dimer) == len(mono2):
    for i in range(0, len(dimer)):
        input_names.append([dimer[i],mono1[i],mono2[i]])
else:
    print("ERROR: Something's wrong with the input ref list!!!")
    exit()

#print input_names

# Call to read parameters file
# line_opt contains these lines need to be optimized
# line_no_opt are the rest

line_no_opt, line_opt, parm_opt_list, parm_position, stype_position, stype_index, n_pound_opt = read_parm(parm_file)


#print parm_opt_list

parm_opt = numpy.zeros(len(parm_opt_list))

for i in range(0, len(parm_opt_list)):
    parm_opt[i] = float(parm_opt_list[i])

#print parm_opt

if idebug == 1 :
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("Debug mode ")
    print("")
    print("The TINKER executable is:     ", tinker_exe)
    print("The energy type to calculate: ", ene_type)
    print("")
    print("Current working dir is:       ", working_dir)
    print("The coordinates dir is:       ", coord_dir)
    print("The input key file is:        ", parm_file)
    print("")
    print("The input list of reference xyz, QM energy and weight is: ", ref_file)
    print("")
    print("The output dir is:            ", out_dir)
    print("The output file is:           ", out_file)
    print("Number of pairs in the data set:                  ", len(input_names))
    print("Number of parallel processes:                     ", nproc)
    print("Number of parameters assigned the #opt:           ", n_pound_opt)
    print("Actual number of parameters need to be optimized: ", len(parm_opt_list))
    print("")
    print("List of parameters need to be optmized: ")
    print("---> ", parm_opt)
    print("")
    print("Parameters bounded?                               ", if_bound)
    if if_bound :
        #print "Boundary ---> ", str(bound_list)
        print("   --------------------------------------------------------")
        print("   || %16s %16s %16s ||" % ("parameter", "lower_bound", "upper_bound"))
        for i in range(0, len(parm_opt_list)) :
            print("   || %16.5f %16.5f %16.5f ||" % (float(parm_opt_list[i]), float(bound_list[i][0]), float(bound_list[i][1])))
        print("   --------------------------------------------------------")
    print("")
    print("The optimizer is :                                      ", optimizer, " : ", optimizer_name[optimizer])
    print("")
    print("The gradient tolerance to converge the optimization is: ", gtol)
    print("")
    print("The step size to evaluate the numerical gradient is:    ", delta_parm)
    print("")
    print("Start the optimization: ")
    print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    print("")

f_out.write( "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n" )
f_out.write( "The TINKER executable is:     " + tinker_exe + "\n" )
f_out.write( "The energy type to calculate: " + ene_type + "\n" )
f_out.write( "" + "\n" )
f_out.write( "Current working dir is:       " + working_dir + "\n" )
f_out.write( "The coordinates dir is:       " + coord_dir + "\n" )
f_out.write( "The input key file is:        " + parm_file + "\n" )
f_out.write( "" + "\n" )
f_out.write( "The input list of reference xyz, QM energy and weight is: " + ref_file + "\n" )
f_out.write( "" + "\n" )
f_out.write( "The output dir is:            " + out_dir + "\n" )
f_out.write( "The output file is:           " + out_file + "\n" )
f_out.write( "" + "\n" )
f_out.write( "Number of pairs in the data set:                  " + str(len(input_names)) + "\n" )
f_out.write( "Number of parallel processes:                     " + str(nproc) + "\n" )
f_out.write( "Number of parameters assigned the #opt:           " + str(n_pound_opt) + "\n" )
f_out.write( "Actual number of parameters need to be optimized: " + str(len(parm_opt_list)) + "\n" )
f_out.write( "" + "\n" )
f_out.write( "List of parameters need to be optmized: " + "\n" )
f_out.write( "---> " + str(parm_opt) + "\n" )
f_out.write( "" + "\n" )
f_out.write( "Parameters bounded?                               " + str(if_bound) + "\n" )
if if_bound :
    #f_out.write( "Boundary ---> " + str(bound_list) + "\n" )
    f_out.write( "   --------------------------------------------------------" + "\n" )
    f_out.write( "   || %16s %16s %16s ||" % ("parameter", "lower_bound", "upper_bound") + "\n" )
    for i in range(0, len(parm_opt_list)):
        f_out.write( "   || %16.5f %16.5f %16.5f ||" % (float(parm_opt_list[i]), float(bound_list[i][0]), float(bound_list[i][1])) + "\n" )
    f_out.write( "   --------------------------------------------------------" + "\n" )
f_out.write( "" + "\n" )
f_out.write( "The optimizer is :                                      " + str(optimizer) + " : " + optimizer_name[optimizer] + "\n" )
f_out.write( "The gradient tolerance to converge the optimization is: " + str(gtol) + "\n" )
f_out.write( "" + "\n" )
f_out.write( "The step size to evaluate the numerical gradient is:    " + str(delta_parm) + "\n" )
f_out.write( "" + "\n" )
f_out.write( "Start the optimization: " + "\n" )
f_out.write( "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" + "\n" )
f_out.write( "" + "\n" )

##########################
xtol = gtol
ftol = xtol

# Call to optimizer
if optimizer == 0 :
    #if idebug == 1 :
    #    print "OUT: Using simplex algorithm first..."
    #    print ""

    #f_out.write( "OUT: Using simplex algorithm first...\n" )
    #f_out.write( "\n" )

    ## use fmin to minimize it first before goes to fmin_bfgs
    #tmp=optimize.fmin(cost_func, parm_opt, (input_names, "", "ene_call"), xtol=0.01, ftol=0.01, maxiter=None, maxfun=len(parm_opt_list)*3, full_output=0)

    ##print tmp
    ##print "ggggggg", parm_opt
    #for i in xrange(0, len(parm_opt)) :
    #    parm_opt[i] = float(tmp[i])

    #if idebug == 1 :
    #    print "OUT: Starting fmin_bfgs with the parm set:"
    #    print "---> ", parm_opt
    #    print ""

    #f_out.write( "OUT: Starting fmin_bfgs with the parm set:\n" )
    #f_out.write( "---> " + str(parm_opt) + "\n" )
    #f_out.write( "\n" )

    # use fmin_bfgs
    my_optimize.fmin_bfgs(cost_func, parm_opt, num_grad, (input_names, "", "ene_call"), gtol)

elif optimizer == 1 :
    # use fmin 
    my_optimize.fmin(cost_func, parm_opt, (input_names, "", "ene_call"), xtol=0.001, ftol=0.001, maxiter=None, maxfun=None, full_output=0)
elif optimizer == 2 :
    # use fmin_cg
    #my_optimize.fmin_cg(cost_func, parm_opt, num_grad, (input_names, "", "ene_call"), gtol)
    print("Not verified yet!!!")
    print("Exit!!!")
    exit()
elif optimizer == 3 :
    # use fmin_powell
    my_optimize.fmin_powell(cost_func, parm_opt, (input_names, "", "ene_call"), xtol=0.001, ftol=0.001, maxiter=None, maxfun=None, full_output=0)
else :
    print("Unknown flag for optimizer !!!")
    print("Exit!!!")
    help_message()
    exit()

seconds = 0.0
minutes = 0.0
hours = 0.0
seconds = time.time() - start_time
if seconds > 60 :
    minutes, seconds = divmod(seconds, 60)
if minutes > 60 :
    hours, minutes = divmod(minutes, 60)

if (idebug == 1) : 
    print("")
    print("Wall time : %10i hours %3i minutes %5.1f seconds" % (int(hours), int(minutes), seconds))
f_out.write( "" + "\n" )
f_out.write( "Wall time : %10i hours %3i minutes %5.1f seconds" % (int(hours), int(minutes), seconds) + "\n" )

#optimize.fmin_bfgs(cost_func, parm_opt, num_grad, (input_names, "", "ene_call"), gtol=0.01)
#result=optimize.fmin_bfgs(cost_func, parm_opt, num_grad, (input_names, ""), gtol=0.3)
#result=optimize.fmin_bfgs(cost_func, parm_opt, None, (input_names, ""), gtol=0.1)
#print result

f_out.close()

