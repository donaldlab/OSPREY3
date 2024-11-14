import numpy as np
import pickle
import pathminus_ph
from pathlib import Path

model_path = Path(__file__).parent / 'path_minus_weights.pkl'

with open(model_path, "rb") as f:
    model = pickle.load(f)

def pathminus_predict(protein_file, ligand_file):
    images = np.array(pathminus_ph.get_all_images(protein_file, ligand_file))
    images_reshaped = np.array(images).reshape(1, images.size)
    prediction = model.predict(images_reshaped)[0]
    return prediction

if __name__ == '__main__':
    import time
    begin_time = time.time()

    protein_file = Path(__file__).parent / '../example/1a1e_protein.pdb'
    ligand_file = Path(__file__).parent / '../example/1a1e_ligand.mol2'
    binding_score = pathminus_predict(protein_file, ligand_file)

    print(f'Binding score: {binding_score:.3f}')
    print(f'Time taken: {time.time() - begin_time:.3f} seconds')

