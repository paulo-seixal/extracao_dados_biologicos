from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.metrics import classification_report, confusion_matrix, matthews_corrcoef, precision_score
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import BernoulliNB, MultinomialNB
from sklearn.linear_model import LogisticRegression
import pandas as pd
import matplotlib.pyplot as plt

def bookmaker(cm):
    tn, fp, fn, tp = cm.ravel()
    return (tp / (tp + fn)) - (fp / (fp + tn))

def fowlkes_mallows(confusion_matrix):
    import numpy as np

    tn, fp, fn, tp = confusion_matrix.ravel()
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    return np.sqrt(precision * recall)




def svm(fingerprints, Xs_train, Xs_test, y_train, test_set_final, train_set_descriptors, test_set_descriptors):
    '''
    Função para aplicação do modelo Linear SVM com otimização de hiperparametros
    '''

    for fp, X_train, X_test in zip(fingerprints,Xs_train, Xs_test):
        param_grid = {'C':[0.1, 1, 10, 100, 1000],
                    'loss':['hinge','squared_hinge'],
                    'dual':[True, False]}
        

        # Classificador SVM
        classifier = GridSearchCV(LinearSVC(), param_grid, cv = 5) 
        classifier.fit(X_train, y_train)

        #Implementar o modelo treinado na previsão da  mutagenicidade
        y_pred = classifier.predict(X_test) 
        y_test = test_set_final['Y'] 
        
        print(f'##### {fp} #####')
        
        # Classification Report
        print(classification_report(y_test, y_pred))

        # Matriz da confusão
        cm = confusion_matrix(y_test, y_pred)
        print("Confusion Matrix ({}) :".format(fp))
        print(cm)
        print()
        print("Cross-Validation Accuracy ({}) : {:.2f}%".format(fp, classifier.cv_results_['mean_test_score'][classifier.best_index_] * 100))

        
        # Cálculo do MCC
        mcc = matthews_corrcoef(y_test, y_pred)
        print("MCC ({}) : {:.2f}".format(fp, mcc))
        
        # Cálculo do Fowlkes-Mallows index
        fm = fowlkes_mallows(cm)
        print("Fowlkes-Mallows Index ({}) : {:.2f}".format(fp, fm))
        
        # Cálculo do Bookmaker informedness
        bm = bookmaker(cm)
        print("Bookmaker Informedness ({}) : {:.2f}".format(fp, bm))
        
        # Cálculo do deltaP (Δp)
        delta_p = precision_score(y_test, y_pred) - sum(y_test) / len(y_test)
        print("DeltaP ({}) : {:.2f}".format(fp, delta_p))
        print('\n##########################################################\n')


    param_grid = {'C':[0.1, 1, 10, 100, 1000],
                    'loss':['hinge','squared_hinge'],
                    'dual':[True, False]}
    scaler = StandardScaler()

    X_train_descriptors = pd.DataFrame(scaler.fit_transform(train_set_descriptors.iloc[:,:-1]))
    y_train_descriptors = train_set_descriptors['Y']
    X_test_descriptors = pd.DataFrame(scaler.fit_transform(test_set_descriptors.iloc[:,:-1]))
    y_test_descriptors = test_set_descriptors['Y']
    # Classificador SVM
    classifier = GridSearchCV(LinearSVC(), param_grid, cv = 5) 
    classifier.fit(X_train_descriptors, y_train_descriptors)

    #Implementar o modelo treinado na previsão da  mutagenicidade
    y_pred_descriptors = classifier.predict(X_test_descriptors) 

    print(f'##### Descriptors #####')

    # Classification Report
    print(classification_report(y_test_descriptors, y_pred_descriptors))

    # Matriz da confusão
    cm = confusion_matrix(y_test_descriptors, y_pred_descriptors)
    print("Confusion Matrix (Descriptors) :")
    print(cm)
    print()
    print("Cross-Validation Accuracy (Descriptors) : {:.2f}%".format(classifier.cv_results_['mean_test_score'][classifier.best_index_] * 100))

    # Cálculo do MCC
    mcc = matthews_corrcoef(y_test, y_pred)
    print("MCC (Descriptors) : {:.2f}".format(mcc))

    # Cálculo do Fowlkes-Mallows index
    fm = fowlkes_mallows(cm)
    print("Fowlkes-Mallows Index (Descriptors) : {:.2f}".format(fm))

    # Cálculo do Bookmaker informedness
    bm = bookmaker(cm)
    print("Bookmaker Informedness (Descriptors) : {:.2f}".format(bm))

    # Cálculo do deltaP (Δp)
    delta_p = precision_score(y_test_descriptors, y_pred_descriptors) - sum(y_test_descriptors) / len(y_test_descriptors)
    print("DeltaP (Descriptors) : {:.2f}".format(delta_p))
    print('\n##########################################################\n')


def svm_curves(fingerprints, train_set_descriptors, test_set_descriptors, Xs_train, Xs_test, y_train, test_set_final):
    '''
    Representação gráfica das curvas ROC e PR do modelo SVM
    '''
    from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay

    scaler = StandardScaler()
    X_train_descriptors = pd.DataFrame(scaler.fit_transform(train_set_descriptors.iloc[:,:-1]))
    y_train_descriptors = train_set_descriptors['Y']
    X_test_descriptors = pd.DataFrame(scaler.fit_transform(test_set_descriptors.iloc[:,:-1]))
    y_test_descriptors = test_set_descriptors['Y'] 
    y_test = test_set_final['Y'] 

    fig, [ax_roc, ax_det] = plt.subplots(1, 2, figsize=(11, 5))

    for fp, X_train, X_test in zip(fingerprints, Xs_train, Xs_test):
        # Classifier SVM
        classifier = LinearSVC()
        classifier.fit(X_train, y_train)

        RocCurveDisplay.from_estimator(classifier, X_test, y_test, ax=ax_roc, name=fp)
        PrecisionRecallDisplay.from_estimator(classifier, X_test, y_test, ax=ax_det, name=fp)
        

    classifier = LinearSVC()
    classifier.fit(X_train_descriptors, y_train_descriptors)

    RocCurveDisplay.from_estimator(classifier, X_test_descriptors, y_test_descriptors, ax=ax_roc, name='Descriptors')
    PrecisionRecallDisplay.from_estimator(classifier, X_test_descriptors, y_test_descriptors, ax=ax_det, name='Descriptors')


    ax_roc.set_title("Curvas ROC - SVM")
    ax_det.set_title("Curvas PR - SVM")

    ax_roc.grid(linestyle="--")
    ax_det.grid(linestyle="--")

    plt.legend()
    plt.show()



def knn(fingerprints, Xs_train, Xs_test, y_train, test_set_final, train_set_descriptors, test_set_descriptors):
    '''
    Função para aplicação do modelo KNN com otimização de hiperparametros
    '''
    param_grid = {'n_neighbors':[2, 5, 10, 20, 30],
                    'weights':['uniform','distance']}

    for fp, X_train, X_test in zip(fingerprints,Xs_train, Xs_test):
               

        # Classificador SVM
        classifier = GridSearchCV(KNeighborsClassifier(), param_grid, cv = 5) 
        classifier.fit(X_train, y_train)

        #Implementar o modelo treinado na previsão da  mutagenicidade
        y_pred = classifier.predict(X_test) 
        y_test = test_set_final['Y'] 
        
        print(f'##### {fp} #####')
        
        # Classification Report
        print(classification_report(y_test, y_pred))

        # Matriz da confusão
        cm = confusion_matrix(y_test, y_pred)
        print("Confusion Matrix ({}) :".format(fp))
        print(cm)
        print()
        print("Cross-Validation Accuracy ({}) : {:.2f}%".format(fp, classifier.cv_results_['mean_test_score'][classifier.best_index_] * 100))

        
        # Cálculo do MCC
        mcc = matthews_corrcoef(y_test, y_pred)
        print("MCC ({}) : {:.2f}".format(fp, mcc))
        
        # Cálculo do Fowlkes-Mallows index
        fm = fowlkes_mallows(cm)
        print("Fowlkes-Mallows Index ({}) : {:.2f}".format(fp, fm))
        
        # Cálculo do Bookmaker informedness
        bm = bookmaker(cm)
        print("Bookmaker Informedness ({}) : {:.2f}".format(fp, bm))
        
        # Cálculo do deltaP (Δp)
        delta_p = precision_score(y_test, y_pred) - sum(y_test) / len(y_test)
        print("DeltaP ({}) : {:.2f}".format(fp, delta_p))
        print('\n##########################################################\n')


    scaler = StandardScaler()

    X_train_descriptors = pd.DataFrame(scaler.fit_transform(train_set_descriptors.iloc[:,:-1]))
    y_train_descriptors = train_set_descriptors['Y']
    X_test_descriptors = pd.DataFrame(scaler.fit_transform(test_set_descriptors.iloc[:,:-1]))
    y_test_descriptors = test_set_descriptors['Y']
    # Classificador SVM
    classifier = GridSearchCV(KNeighborsClassifier(), param_grid, cv = 5) 
    classifier.fit(X_train_descriptors, y_train_descriptors)

    #Implementar o modelo treinado na previsão da  mutagenicidade
    y_pred_descriptors = classifier.predict(X_test_descriptors) 

    print('##### Descriptors #####')

    # Classification Report
    print(classification_report(y_test_descriptors, y_pred_descriptors))

    # Matriz da confusão
    cm = confusion_matrix(y_test_descriptors, y_pred_descriptors)
    print("Confusion Matrix (Descriptors) :")
    print(cm)
    print()
    print("Cross-Validation Accuracy (Descriptors) : {:.2f}%".format(classifier.cv_results_['mean_test_score'][classifier.best_index_] * 100))

    # Cálculo do MCC
    mcc = matthews_corrcoef(y_test, y_pred)
    print("MCC (Descriptors) : {:.2f}".format(mcc))

    # Cálculo do Fowlkes-Mallows index
    fm = fowlkes_mallows(cm)
    print("Fowlkes-Mallows Index (Descriptors) : {:.2f}".format(fm))

    # Cálculo do Bookmaker informedness
    bm = bookmaker(cm)
    print("Bookmaker Informedness (Descriptors) : {:.2f}".format(bm))

    # Cálculo do deltaP (Δp)
    delta_p = precision_score(y_test_descriptors, y_pred_descriptors) - sum(y_test_descriptors) / len(y_test_descriptors)
    print("DeltaP (Descriptors) : {:.2f}".format(delta_p))
    print('\n##########################################################\n')


def knn_curves(fingerprints, train_set_descriptors, test_set_descriptors, Xs_train, Xs_test, y_train, test_set_final):
    '''
    Representação gráfica das curvas ROC e PR do modelo KNN
    '''
    from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay

    scaler = StandardScaler()
    X_train_descriptors = pd.DataFrame(scaler.fit_transform(train_set_descriptors.iloc[:,:-1]))
    y_train_descriptors = train_set_descriptors['Y']
    X_test_descriptors = pd.DataFrame(scaler.fit_transform(test_set_descriptors.iloc[:,:-1]))
    y_test_descriptors = test_set_descriptors['Y'] 
    y_test = test_set_final['Y'] 

    fig, [ax_roc, ax_det] = plt.subplots(1, 2, figsize=(11, 5))

    for fp, X_train, X_test in zip(fingerprints, Xs_train, Xs_test):
        # Classifier knn
        classifier = KNeighborsClassifier()
        classifier.fit(X_train, y_train)

        RocCurveDisplay.from_estimator(classifier, X_test, y_test, ax=ax_roc, name=fp)
        PrecisionRecallDisplay.from_estimator(classifier, X_test, y_test, ax=ax_det, name=fp)
        

    classifier = KNeighborsClassifier()
    classifier.fit(X_train_descriptors, y_train_descriptors)

    RocCurveDisplay.from_estimator(classifier, X_test_descriptors, y_test_descriptors, ax=ax_roc, name='Descriptors')
    PrecisionRecallDisplay.from_estimator(classifier, X_test_descriptors, y_test_descriptors, ax=ax_det, name='Descriptors')


    ax_roc.set_title("Curvas ROC - KNN")
    ax_det.set_title("Curvas PR - KNN")

    ax_roc.grid(linestyle="--")
    ax_det.grid(linestyle="--")

    plt.legend()
    plt.show()



def naive_bayes(fingerprints, Xs_train, Xs_test, y_train, test_set_final, train_set_descriptors, test_set_descriptors):
    '''
    Função para aplicação do modelo Naive Bayes com otimização de hiperparametros
    '''
    param_grid = {'alpha': [0, 0.01, 0.5, 1.0, 5, 10],
              'force_alpha':[True, False],
              'fit_prior': [True, False]}

    for fp, X_train, X_test in zip(fingerprints,Xs_train, Xs_test):
               

        # Classificador SVM
        classifier = GridSearchCV(BernoulliNB(), param_grid, cv = 5) 
        classifier.fit(X_train, y_train)

        #Implementar o modelo treinado na previsão da  mutagenicidade
        y_pred = classifier.predict(X_test) 
        y_test = test_set_final['Y'] 
        
        print(f'##### {fp} #####')
        
        # Classification Report
        print(classification_report(y_test, y_pred))

        # Matriz da confusão
        cm = confusion_matrix(y_test, y_pred)
        print("Confusion Matrix ({}) :".format(fp))
        print(cm)
        print()
        print("Cross-Validation Accuracy ({}) : {:.2f}%".format(fp, classifier.cv_results_['mean_test_score'][classifier.best_index_] * 100))

        
        # Cálculo do MCC
        mcc = matthews_corrcoef(y_test, y_pred)
        print("MCC ({}) : {:.2f}".format(fp, mcc))
        
        # Cálculo do Fowlkes-Mallows index
        fm = fowlkes_mallows(cm)
        print("Fowlkes-Mallows Index ({}) : {:.2f}".format(fp, fm))
        
        # Cálculo do Bookmaker informedness
        bm = bookmaker(cm)
        print("Bookmaker Informedness ({}) : {:.2f}".format(fp, bm))
        
        # Cálculo do deltaP (Δp)
        delta_p = precision_score(y_test, y_pred) - sum(y_test) / len(y_test)
        print("DeltaP ({}) : {:.2f}".format(fp, delta_p))
        print('\n##########################################################\n')


    scaler = StandardScaler()

    X_train_descriptors = pd.DataFrame(scaler.fit_transform(train_set_descriptors.iloc[:,:-1]))
    y_train_descriptors = train_set_descriptors['Y']
    X_test_descriptors = pd.DataFrame(scaler.fit_transform(test_set_descriptors.iloc[:,:-1]))
    y_test_descriptors = test_set_descriptors['Y']
    # Classificador SVM
    classifier = GridSearchCV(BernoulliNB(), param_grid, cv = 5) 
    classifier.fit(X_train_descriptors, y_train_descriptors)

    #Implementar o modelo treinado na previsão da  mutagenicidade
    y_pred_descriptors = classifier.predict(X_test_descriptors) 

    print('##### Descriptors #####')

    # Classification Report
    print(classification_report(y_test_descriptors, y_pred_descriptors))

    # Matriz da confusão
    cm = confusion_matrix(y_test_descriptors, y_pred_descriptors)
    print("Confusion Matrix (Descriptors) :")
    print(cm)
    print()
    print("Cross-Validation Accuracy (Descriptors) : {:.2f}%".format(classifier.cv_results_['mean_test_score'][classifier.best_index_] * 100))

    # Cálculo do MCC
    mcc = matthews_corrcoef(y_test, y_pred)
    print("MCC (Descriptors) : {:.2f}".format(mcc))

    # Cálculo do Fowlkes-Mallows index
    fm = fowlkes_mallows(cm)
    print("Fowlkes-Mallows Index (Descriptors) : {:.2f}".format(fm))

    # Cálculo do Bookmaker informedness
    bm = bookmaker(cm)
    print("Bookmaker Informedness (Descriptors) : {:.2f}".format(bm))

    # Cálculo do deltaP (Δp)
    delta_p = precision_score(y_test_descriptors, y_pred_descriptors) - sum(y_test_descriptors) / len(y_test_descriptors)
    print("DeltaP (Descriptors) : {:.2f}".format(delta_p))
    print('\n##########################################################\n')


def naive_bayes_curves(fingerprints, train_set_descriptors, test_set_descriptors, Xs_train, Xs_test, y_train, test_set_final):
    '''
    Representação gráfica das curvas ROC e PR do modelo Naive Bayes
    '''
    from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay

    scaler = StandardScaler()
    X_train_descriptors = pd.DataFrame(scaler.fit_transform(train_set_descriptors.iloc[:,:-1]))
    y_train_descriptors = train_set_descriptors['Y']
    X_test_descriptors = pd.DataFrame(scaler.fit_transform(test_set_descriptors.iloc[:,:-1]))
    y_test_descriptors = test_set_descriptors['Y'] 
    y_test = test_set_final['Y'] 

    fig, [ax_roc, ax_det] = plt.subplots(1, 2, figsize=(11, 5))

    for fp, X_train, X_test in zip(fingerprints, Xs_train, Xs_test):
        # Classifier knn
        classifier = BernoulliNB()
        classifier.fit(X_train, y_train)

        RocCurveDisplay.from_estimator(classifier, X_test, y_test, ax=ax_roc, name=fp)
        PrecisionRecallDisplay.from_estimator(classifier, X_test, y_test, ax=ax_det, name=fp)
        

    classifier = BernoulliNB()
    classifier.fit(X_train_descriptors, y_train_descriptors)

    RocCurveDisplay.from_estimator(classifier, X_test_descriptors, y_test_descriptors, ax=ax_roc, name='Descriptors')
    PrecisionRecallDisplay.from_estimator(classifier, X_test_descriptors, y_test_descriptors, ax=ax_det, name='Descriptors')


    ax_roc.set_title("Curvas ROC - Naive Bayes")
    ax_det.set_title("Curvas PR - Naive Bayes")

    ax_roc.grid(linestyle="--")
    ax_det.grid(linestyle="--")

    plt.legend()
    plt.show()

def dl_model_metrics(model, X_test, y_test, title):
    '''
    Função para avaliação dos modelos deep learning, através de várias métricas
    '''
    y_pred = model.predict(X_test)

    loss, accuracy = model.evaluate(X_test, y_test)

    # Apply a threshold to convert probabilities to binary values
    threshold = 0.5
    y_pred_binary = (y_pred >= threshold).astype(int)
    y_test_binary = (y_test >= threshold).astype(int)

    print(f'\n################## {title} ##################\n')
    print('Test loss:', loss)
    print('Test accuracy:', accuracy)

    print('\n- - - - - - - - - - - - - - - - - - - - - - - - -')
    # Confusion matrix
    cm = confusion_matrix(y_test_binary, y_pred_binary)
    print("\nConfusion Matrix:")
    print(cm)
    print('- - - - - - - - - - - - - - - - - - - - - - - - -\n')

    # classification report
    print("Classification Report:\n", classification_report(y_test_binary, y_pred_binary))
    print('- - - - - - - - - - - - - - - - - - - - - - - - -\n')
    # Calculate MCC
    mcc = matthews_corrcoef(y_test_binary, y_pred_binary)
    print("MCC: {:.2f}".format(mcc))

    # Calculate Fowlkes-Mallows index
    fm = fowlkes_mallows(cm)
    print("Fowlkes-Mallows Index: {:.2f}".format(fm))

    # Cálculo do Bookmaker informedness
    bm = bookmaker(cm)
    print("Bookmaker Informedness: {:.2f}".format(bm))

    # Cálculo do deltaP (Δp)
    delta_p = precision_score(y_test_binary, y_pred_binary) - sum(y_test) / len(y_test_binary)
    print("DeltaP: {:.2f}".format(delta_p))


def plot_acc_loss_dl(model, title):
    '''
    Função para representação gráfica dos valores de loss e accuracy no treino e na validação dos modelos deep learning
    '''
    # Obter as métricas de treino
    train_loss = model.history['loss']
    train_acc = model.history['accuracy']

    # Obter as métricas de validação
    val_loss = model.history['val_loss']
    val_acc = model.history['val_accuracy']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

    # Criar o gráfico
    epochs = range(1, len(train_loss) + 1)
    ax1.plot(epochs, train_loss, 'b', label='Training Loss')
    ax1.plot(epochs, val_loss, 'r', label='Validation Loss')
    ax1.set_title(f'Training and Validation Loss - {title}')
    ax1.set_xlabel('Epochs')
    ax1.set_ylabel('Loss')
    ax1.legend()

    # Criar o gráfico de acurácia
    ax2.plot(epochs, train_acc, 'b', label='Training Accuracy')
    ax2.plot(epochs, val_acc, 'r', label='Validation Accuracy')
    ax2.set_title(f'Training and Validation Accuracy - {title}')
    ax2.set_xlabel('Epochs')
    ax2.set_ylabel('Accuracy')
    ax2.legend()
