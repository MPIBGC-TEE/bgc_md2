from os import environ
from datetime import datetime

from dask_mpi import initialize
from distributed import Client

from sklearn import datasets
from sklearn.model_selection import train_test_split

# Get the Dask version of GridSearchCV
from dask_ml.model_selection import GridSearchCV

from sklearn.metrics import classification_report
from sklearn.svm import SVC


def run_test(client):
    print(__doc__)

    # Loading the Digits dataset
    digits = datasets.load_digits()

    # To apply an classifier on this data, we need to flatten the image, to
    # turn the data in a (samples, feature) matrix:
    n_samples = len(digits.images)
    X = digits.images.reshape((n_samples, -1))
    y = digits.target

    # Split the dataset in two equal parts
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.5, random_state=0)

    # Set the parameters by cross-validation
    tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
                         'C': [1, 10, 100, 1000]},
                        {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]

    scores = ['precision', 'recall']

    for score in scores:
        print("# Tuning hyper-parameters for %s" % score)
        print()

        # scheduler=client makes sure that Dask uses the correct communications
        clf = GridSearchCV(
            SVC(), tuned_parameters, scoring='%s_macro' % score,
            scheduler=client
        )
        clf.fit(X_train, y_train)

        print("Best parameters set found on development set:")
        print()
        print(clf.best_params_)
        print()
        print("Grid scores on development set:")
        print()
        means = clf.cv_results_['mean_test_score']
        stds = clf.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, clf.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))
        print()

        print("Detailed classification report:")
        print()
        print("The model is trained on the full development set.")
        print("The scores are computed on the full evaluation set.")
        print()
        y_true, y_pred = y_test, clf.predict(X_test)
        print(classification_report(y_true, y_pred))
        print()

    # Note the problem is too easy: the hyperparameter plateau is too flat and
    # the output model is the same for precision and recall with ties in
    # quality.


def main():
    # Work out from the environment how many threads to allocate
    num_threads = int(environ.get(
        'SLURM_CPUS_PER_TASK',
        environ.get('OMP_NUM_THREADS', 1)
    ))

    # Create the Dask workers
    initialize(interface='ib0', nthreads=num_threads)

    # Create the Dask object that will manage the communications
    client = Client()

    start = datetime.now()
    run_test(client=client)
    end = datetime.now()

    print(f"Time taken: {end - start}")


if __name__ == '__main__':
    main()