def get_v_cholerae_oric():
    with open('data/v_cholerae_oric.txt', 'r') as file:
        content = file.read()
        return content
    
def get_v_cholerae():
    with open('data/Vibrio_cholerae.txt', 'r') as file:
        content = file.read()
        return content
    
def get_t_petrophila_oric():
    with open('data/t_petrophila_oriC.txt', 'r') as file:
        content = file.read()
        return content
    
def get_t_petrophila():
    with open('data/Thermotoga_petrophila.txt', 'r') as file:
        content = file.read()
        return content
    
def get_e_coli():
    with open('data/E_coli.txt', 'r') as file:
        content = file.read()
        return content
    
def get_dosr_motif():
    with open('data/DosR.txt', 'r') as file:
        content = file.read()
        return content.split('\n')